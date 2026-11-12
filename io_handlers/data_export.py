from pathlib import Path
import json
import pandas as pd
from core.models import AnalysisResult


def export_result(result: AnalysisResult, output_dir: str | Path, min_score_threshold: float = 0.12) -> dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    
    json_path = out / f"{result.case_id}_analysis.json"
    json_path.write_text(result.model_dump_json(indent=2), encoding="utf-8")
    
    tsv_path = out / f"{result.case_id}_predictions.tsv"
    
    best_prediction = "no_causal_variant_in_vcf"
    
    if result.candidates:
        top_candidate = result.candidates[0]
        
        
        if top_candidate.score >= min_score_threshold:
            best_variant_id = None
            
            
            if getattr(top_candidate, 'variant_id', None):
                best_variant_id = top_candidate.variant_id
            else:
                
                gene_variants = [v for v in result.variants if v.gene == top_candidate.gene]
                if gene_variants:
                    
                    best_variant_id = gene_variants[0].variant_id
            
            if best_variant_id:
                best_prediction = best_variant_id
                
    rows = [{
        "case_id": result.case_id,
        "prediction": best_prediction
    }]
    
    df = pd.DataFrame(rows, columns=["case_id", "prediction"])
    df.to_csv(tsv_path, sep="\t", index=False)
    
    return {"json": str(json_path), "predictions_tsv": str(tsv_path)}