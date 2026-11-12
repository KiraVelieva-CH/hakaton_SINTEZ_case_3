from pathlib import Path
import pandas as pd
from .base_analyzer import BaseAnalyzer

class PhenotypeAnalyzer(BaseAnalyzer):
    name = "phenotype"
    def __init__(self, resource_path=None):
        self.associations = {}
        if resource_path and Path(resource_path).exists():
            df = pd.read_csv(resource_path, sep="\t").fillna("")
            # Expected columns: gene, hpo_id (optional term), evidence (optional)
            for _, r in df.iterrows():
                gene = str(r.get("gene", "")).strip()
                key = str(r.get("hpo_id", r.get("term", ""))).strip().lower()
                if gene and key: self.associations.setdefault(gene, set()).add(key)

    def analyze(self, context):
        phenos = context["phenotypes"]
        keys = {p.hpo_id.lower() for p in phenos if p.hpo_id}
        keys |= {p.term.lower() for p in phenos}
        genes = {v.gene for v in context["variants"] if v.gene}
        results = {}
        for gene in genes:
            assoc = self.associations.get(gene, set())
            hits = sorted(keys & assoc)
            score = len(hits) / max(1, len(keys)) if assoc else 0.0
            results[gene] = {"score": score, "matched_phenotypes": hits,
                             "known_associations": len(assoc), "method": "local_gene_phenotype_overlap"}
        return {"genes": results, "note": "No local association file means phenotype scores are zero, not evidence against a gene."}
