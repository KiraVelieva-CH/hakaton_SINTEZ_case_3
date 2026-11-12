from .base_analyzer import BaseAnalyzer

class InheritanceAnalyzer(BaseAnalyzer):
    name = "inheritance"

    def __init__(self):
        
        self.inheritance_map = {
            "FBN1": "AD", "FGFR2": "AD", "MAP2K1": "AD", 
            "SMAD3": "AD", "HERC1": "AR", "CHST3": "AR", "BBS4": "AR"
        }

    def analyze(self, context):
        genes = {}
        for v in context["variants"]:
            if not v.gene: continue
            
            gt = (v.genotype or "").replace("|", "/")
            alt_count = gt.count("1")
            zygosity = "Гомозигота" if alt_count == 2 else "Гетерозигота" if alt_count == 1 else "Референс"
            
            mode = self.inheritance_map.get(v.gene, "UNKNOWN")
            score = 0.5 
            
            if mode == "AD" and alt_count >= 1:
                score = 1.0
            elif mode == "AR" and alt_count == 2:
                score = 1.0
            elif mode == "AR" and alt_count == 1:
                score = 0.6 
                
            genes.setdefault(v.gene, []).append({
                "score": score,
                "variant_id": v.variant_id,
                "zygosity": zygosity,
                "mode": mode
            })
            
        
        return {"genes": {g: max(vars, key=lambda x: x["score"]) for g, vars in genes.items()}}
