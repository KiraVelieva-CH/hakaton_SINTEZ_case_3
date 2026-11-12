from pathlib import Path
from cyvcf2 import VCF
from .base_analyzer import BaseAnalyzer

class PathogenicityAnalyzer(BaseAnalyzer):
    name = "pathogenicity"

    def __init__(self, clinvar_vcf_path=None):
        self.clinvar = {}
        
        if clinvar_vcf_path and Path(clinvar_vcf_path).exists():
            print("Загрузка ClinVar... (это может занять 1-2 минуты при первом запуске)")
            try:
                
                vcf = VCF(str(clinvar_vcf_path))
                
                
                for variant in vcf:
                    
                    clnsig = variant.INFO.get('CLNSIG')
                    if not clnsig:
                        continue
                    
                    chrom = variant.CHROM if variant.CHROM.startswith('chr') else f"chr{variant.CHROM}"
                    
                    for alt in variant.ALT:
                        var_id = f"{chrom}:{variant.POS}:{variant.REF}:{alt}"
                        
                        
                        sig = str(clnsig).lower()
                        if "pathogenic" in sig and "conflicting" not in sig:
                            score = 1.0
                            text = "Pathogenic"
                        elif "likely_pathogenic" in sig:
                            score = 0.85
                            text = "Likely pathogenic"
                        elif "benign" in sig:
                            score = 0.0
                            text = "Benign"
                        elif "uncertain" in sig:
                            score = 0.35
                            text = "Uncertain significance"
                        else:
                            continue 
                            
                        self.clinvar[var_id] = {"score": score, "significance": text}
                        
                vcf.close()
                print(f"ClinVar успешно загружен. Найдено {len(self.clinvar)} аннотированных вариантов.")
                
            except Exception as e:
                print(f"Ошибка при чтении ClinVar: {e}")

    def analyze(self, context):
        result = {}
        for v in context["variants"]:
            
            cv_data = self.clinvar.get(v.variant_id)
            
            if cv_data:
                score = cv_data["score"]
                significance = cv_data["significance"]
            else:
                score = None
                significance = "Not found in ClinVar"
                
            
            extras = [x for x in (v.scores.get("alphamissense"), v.scores.get("spliceai")) if x is not None]
            if extras: 
                score = sum(extras) / len(extras) if score is None else 0.7 * score + 0.3 * (sum(extras)/len(extras))
                
            result[v.variant_id] = {
                "score": score, 
                "clinvar": {"significance": significance},
                "functional_scores": v.scores
            }
            
        return {"variants": result, "note": "ClinVar loaded locally."}
