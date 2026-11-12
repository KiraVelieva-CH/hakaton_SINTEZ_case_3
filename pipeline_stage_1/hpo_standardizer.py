from pathlib import Path
import pandas as pd
from core.models import Phenotype

class HPOStandardizer:
    def __init__(self, mapping_path=None):
        self.mapping = {}
        if mapping_path and Path(mapping_path).exists():
            df = pd.read_csv(mapping_path, sep="\t").fillna("")
            for _, r in df.iterrows():
                term = str(r.get("term", "")).strip().lower()
                hid = str(r.get("hpo_id", "")).strip()
                if term and hid: self.mapping[term] = hid

    def standardize(self, phenotypes: list[Phenotype]) -> list[Phenotype]:
        for p in phenotypes:
            if not p.hpo_id:
                p.hpo_id = self.mapping.get(p.term.lower())
        return phenotypes
