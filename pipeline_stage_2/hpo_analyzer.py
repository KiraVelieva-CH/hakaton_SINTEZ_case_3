import csv
from pathlib import Path
from .base_analyzer import BaseAnalyzer

class HPOAnalyzer(BaseAnalyzer):
    name = "phenotype"

    def __init__(self, hpo_mapping_path=None):
        self.gene_to_hpo = {}
        
        if hpo_mapping_path and Path(hpo_mapping_path).exists():
            with open(hpo_mapping_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    gene = row.get("gene_symbol")
                    hpo_id = row.get("hpo_id")
                    if gene and hpo_id:
                        self.gene_to_hpo.setdefault(gene, set()).add(hpo_id)

    def analyze(self, context):
        
        patient_hpos = {p.hpo_id for p in context["phenotypes"] if p.present and p.hpo_id}
        genes = {v.gene for v in context["variants"] if v.gene}
        
        results = {}
        for gene in genes:
            gene_hpos = self.gene_to_hpo.get(gene, set())
            if not gene_hpos or not patient_hpos:
                results[gene] = {"score": 0.0, "matched": 0}
                continue
            
            
            intersection = patient_hpos.intersection(gene_hpos)
            
            
            score = len(intersection) / len(patient_hpos)
            
            results[gene] = {"score": score, "matched": len(intersection)}
            
        return {"genes": results}