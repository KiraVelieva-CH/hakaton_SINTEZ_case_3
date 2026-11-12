from core.models import Gene

def build_gene_ranking(candidates, variants):
    by_gene = {}
    for v in variants:
        if v.gene: by_gene.setdefault(v.gene, []).append(v.variant_id)
    return [Gene(gene_id=c.gene, symbol=c.gene, score=c.score,
                 evidence=[e.statement for e in c.evidence],
                 variants=by_gene.get(c.gene, [])) for c in candidates]
