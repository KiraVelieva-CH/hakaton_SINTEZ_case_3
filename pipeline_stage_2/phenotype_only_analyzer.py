import re
from .base_analyzer import BaseAnalyzer

_HPO_RE = re.compile(r"^HP:\d+$")


class PhenotypeOnlyAnalyzer(BaseAnalyzer):
    """Режим без вариантов: ищем по всему справочнику гены, пересекающиеся с фенотипом."""
    name = "phenotype"

    def __init__(self, genes_to_phenotype_path):
        self.gene_terms = {}
        with open(genes_to_phenotype_path, encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                hpo = next((x.strip().upper() for x in parts if _HPO_RE.match(x.strip().upper())), None)
                gene = parts[1].strip()
                if hpo and gene:
                    self.gene_terms.setdefault(gene, set()).add(hpo)

    def analyze(self, context):
        patient = {p.hpo_id.strip().upper() for p in context["phenotypes"] if p.hpo_id and p.present}
        if not patient:
            return {"genes": {}, "note": "Нет HPO-идентификаторов у присутствующих признаков."}
        genes = {}
        for gene, terms in self.gene_terms.items():
            hits = patient & terms
            if hits:
                genes[gene] = {"score": len(hits) / len(patient), "matched": len(hits),
                               "matched_phenotypes": sorted(hits),
                               "method": "phenotype_only_overlap"}
        return {"genes": genes,
                "note": "Режим без вариантов: результат это список генов, конкретный вариант не предсказывается."}