from core.models import Candidate, EvidenceItem

def rank_candidates(context, block_results, config):
    weights = config.get("weights", {})
    variants_by_gene = {}
    for v in context["variants"]:
        if v.gene: variants_by_gene.setdefault(v.gene, []).append(v)
    
    genes = set(variants_by_gene)
    for block in ("phenotype", "inheritance"):
        genes |= set(block_results.get(block, {}).get("genes", {}))
        
    candidates = []
    for gene in genes:
        
        ph_data = block_results.get("phenotype", {}).get("genes", {}).get(gene, {})
        ph = ph_data.get("score")
        ph_matched = ph_data.get("matched", 0)

        inh_data = block_results.get("inheritance", {}).get("genes", {}).get(gene, {})
        inh = inh_data.get("score")
        zygosity = inh_data.get("zygosity", "неизвестно")
        mode = inh_data.get("mode", "неизвестно")
        
        
        exo = None 
        
        gene_vars = variants_by_gene.get(gene, [])
        
        
        best_variant_id = gene_vars[0].variant_id if gene_vars else None
        path_scores = []
        for v in gene_vars:
            p_score = block_results.get("pathogenicity", {}).get("variants", {}).get(v.variant_id, {}).get("score")
            if p_score is not None:
                path_scores.append(p_score)
        
        path = max(path_scores) if path_scores else None
        
        
        components = {"phenotype": ph, "pathogenicity": path, "inheritance": inh, "exomiser": exo}
        weighted, total_w = 0.0, 0.0
        for k, val in components.items():
            w = float(weights.get(k, 0))
            if val is not None and w > 0:
                weighted += max(0.0, min(1.0, float(val))) * w
                total_w += w
        score = weighted / total_w if total_w else 0.0
        
        
        evidence, caveats = [], []
        
        if ph is not None: 
            evidence.append(EvidenceItem(
                source="Фенотип", 
                statement=f"Совпадение по HPO: {ph:.3f} (найдено {ph_matched} общих симптомов)"
            ))
            
        if inh is not None: 
            evidence.append(EvidenceItem(
                source="Наследование", 
                statement=f"Режим гена: {mode}. Генотип: {zygosity}. Скор: {inh:.3f}"
            ))
            
        if path is not None: 
            evidence.append(EvidenceItem(source="Патогенность", statement=f"Функциональный скор: {path:.3f}"))
            
        if not any(x is not None for x in components.values()): 
            caveats.append("Нет данных для оценки.")
            
        if any(v.filters for v in gene_vars): 
            caveats.append("Есть варианты с флагами QC (низкое качество).")
            
        if ph is None: 
            caveats.append("Нет фенотипической ассоциации.")
            
        if path is None: 
            caveats.append("Нет аннотации патогенности (NA в таблице).")
            
        
        candidates.append(Candidate(
            gene=gene, 
            variant_id=best_variant_id, 
            score=score, 
            component_scores=components,
            evidence=evidence, 
            caveats=caveats
        ))
        
    candidates.sort(key=lambda c: c.score, reverse=True)
    return candidates
