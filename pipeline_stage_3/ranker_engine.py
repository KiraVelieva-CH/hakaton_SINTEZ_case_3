from core.models import Candidate, EvidenceItem

def rank_candidates(context, block_results, config):
    weights = config.get("weights", {})
    variants_by_gene = {}
    for v in context["variants"]:
        if v.gene: variants_by_gene.setdefault(v.gene, []).append(v)
    genes = set(variants_by_gene)
    for block in ("phenotype", "inheritance", "exomiser"):
        genes |= set(block_results.get(block, {}).get("genes", {}))
    candidates = []
    for gene in genes:
        ph = block_results.get("phenotype", {}).get("genes", {}).get(gene, {}).get("score")
        inh = block_results.get("inheritance", {}).get("genes", {}).get(gene, {}).get("score")
        exo = block_results.get("exomiser", {}).get("genes", {}).get(gene, {}).get("score")
        gene_vars = variants_by_gene.get(gene, [])
        path_scores = [block_results.get("pathogenicity", {}).get("variants", {}).get(v.variant_id, {}).get("score") for v in gene_vars]
        path_scores = [x for x in path_scores if x is not None]
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
        if ph is not None: evidence.append(EvidenceItem(source="phenotype", statement=f"Phenotype similarity score: {ph:.3f}"))
        if path is not None: evidence.append(EvidenceItem(source="pathogenicity", statement=f"Best available variant pathogenicity score: {path:.3f}"))
        if inh is not None: evidence.append(EvidenceItem(source="inheritance", statement=f"Genotype heuristic score: {inh:.3f}"))
        if exo is not None: evidence.append(EvidenceItem(source="exomiser", statement=f"Exomiser score: {exo:.3f}"))
        if not any(x is not None for x in components.values()): caveats.append("No scored evidence available.")
        if any(v.filters for v in gene_vars): caveats.append("One or more variants have QC flags; inspect variant-level details.")
        if ph is None: caveats.append("No phenotype association score was available.")
        if path is None: caveats.append("No pathogenicity annotation was available; missing does not mean benign.")
        if inh is not None: caveats.append("Inheritance score is a coarse genotype heuristic, not a clinical inheritance assessment.")
        candidates.append(Candidate(gene=gene, score=score, component_scores=components,
            evidence=evidence, caveats=caveats))
    candidates.sort(key=lambda c: c.score, reverse=True)
    return candidates
