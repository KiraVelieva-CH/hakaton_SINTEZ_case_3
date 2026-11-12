from .base_analyzer import BaseAnalyzer

class InheritanceAnalyzer(BaseAnalyzer):
    name = "inheritance"

    def __init__(self, gene_disease_tsv=None):
        
        self.inheritance_map = {"FBN1": "AD"} 
        
    def analyze(self, context):
        genes = {}
        for v in context["variants"]:
            if not v.gene: continue
            
            
            gt = (v.genotype or "").replace("|", "/")
            alt_count = gt.count("1")
            
            mode = self.inheritance_map.get(v.gene, "UNKNOWN")
            score = 0.0
            
            if mode == "AD":
                
                if alt_count >= 1: score = 1.0
            elif mode == "AR":
                
                if alt_count == 2: score = 1.0
                elif alt_count == 1: score = 0.5 
            else:
                score = 0.5 
                
            genes.setdefault(v.gene, []).append(score)
            
        return {"genes": {g: {"score": max(scores)} for g, scores in genes.items()}}