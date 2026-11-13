from pathlib import Path
from .base_analyzer import BaseAnalyzer

class PathogenicityAnalyzer(BaseAnalyzer):
    name = "pathogenicity"

    def __init__(self, clinvar_vcf_path=None):
        self.clinvar = {}
        self.is_active = False 
        

        if clinvar_vcf_path and Path(clinvar_vcf_path).exists():
            print("Загрузка ClinVar... (это может занять 1-2 минуты при первом запуске)")
            try:

                from cyvcf2 import VCF
                
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
                self.is_active = True
                print(f"ClinVar успешно загружен. Найдено {len(self.clinvar)} аннотированных вариантов.")
                
            except ImportError:
                print("ВНИМАНИЕ: Библиотека cyvcf2 не установлена. Анализ патогенности по ClinVar будет пропущен.")
            except Exception as e:
                print(f"Ошибка при чтении ClinVar: {e}. Анализ патогенности будет пропущен.")
        else:

            print("ClinVar не подключен (файл не указан или отсутствует). Будут использоваться только фенотип и наследование, а также информацию о патогенезе из ваших файлов.")

    def analyze(self, context):
        result = {}
        for v in context["variants"]:
            
            cv_data = self.clinvar.get(v.variant_id)
            
            if cv_data:
                score = cv_data["score"]
                significance = cv_data["significance"]
            else:
                score = None
                significance = "Not found in ClinVar" if self.is_active else "ClinVar disabled"
                

            extras = [x for x in (v.scores.get("alphamissense"), v.scores.get("spliceai")) if x is not None]
            if extras: 
                score = sum(extras) / len(extras) if score is None else 0.7 * score + 0.3 * (sum(extras)/len(extras))
                
            result[v.variant_id] = {
                "score": score, 
                "clinvar": {"significance": significance},
                "functional_scores": v.scores
            }
            

        note = "ClinVar loaded locally." if self.is_active else "ClinVar is disabled. Missing pathogenicity data means UNKNOWN, not benign."
        return {"variants": result, "note": note}