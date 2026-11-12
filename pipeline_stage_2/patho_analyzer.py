from pathlib import Path
import pandas as pd
from .base_analyzer import BaseAnalyzer

class PathogenicityAnalyzer(BaseAnalyzer):
    name = "pathogenicity"
    def __init__(self, clinvar_path=None):
        self.clinvar = {}
        if clinvar_path and Path(clinvar_path).exists():
            df = pd.read_csv(clinvar_path, sep="\t").fillna("")
            for _, r in df.iterrows():
                key = str(r.get("variant_id", "")).strip()
                if key: self.clinvar[key] = {"significance": str(r.get("significance","")), "review_status": str(r.get("review_status",""))}

    def analyze(self, context):
        result = {}
        for v in context["variants"]:
            cv = self.clinvar.get(v.variant_id, {})
            sig = (v.clinvar_significance or cv.get("significance", "")).lower()
            score = None
            if "pathogenic" in sig and "conflict" not in sig: score = 1.0
            elif "likely pathogenic" in sig: score = 0.85
            elif "benign" in sig: score = 0.0
            elif "uncertain" in sig or "vus" in sig: score = 0.35
            # Functional scores are only used when actually supplied in the input table.
            extras = [x for x in (v.scores.get("alphamissense"), v.scores.get("spliceai")) if x is not None]
            if extras: score = sum(extras) / len(extras) if score is None else 0.7 * score + 0.3 * (sum(extras)/len(extras))
            result[v.variant_id] = {"score": score, "clinvar": cv or {"significance": v.clinvar_significance},
                                    "functional_scores": v.scores}
        return {"variants": result, "note": "Missing annotation is unknown; it is not benign. AlphaMissense/SpliceAI are not queried automatically."}
