from __future__ import annotations
from pathlib import Path

from io_handlers.data_loader import load_phenotypes, load_variant_table, load_vcf
from pipeline_stage_2.hpo_analyzer import HPOAnalyzer
from pipeline_stage_2.patho_analyzer import PathogenicityAnalyzer
from pipeline_stage_2.inherit_analyzer import InheritanceAnalyzer
from pipeline_stage_2.phenotype_only_analyzer import PhenotypeOnlyAnalyzer

RESOURCE_KEYS = ("hpo_terms_tsv", "genes_to_phenotype_txt", "clinvar_vcf")

STATUS_RU = {
    "found": "найден", "readable": "прочитан", "missing": "ОТСУТСТВУЕТ",
    "unreadable": "НЕ ПРОЧИТАН", "not_provided": "не задан",
    "ok": "выполнен", "failed": "ОШИБКА", "skipped": "пропущен",
}


def _file_state(path):
    if not path:
        return {"status": "not_provided", "path": None}
    p = Path(path)
    if not p.is_file():
        return {"status": "missing", "path": str(p)}
    return {"status": "found", "path": str(p)}


def _try_load(path, loader):
    state = _file_state(path)
    if state["status"] != "found":
        return [], state
    try:
        data = loader(path)
    except Exception as e:
        state["status"] = "unreadable"
        state["error"] = str(e)
        return [], state
    state["status"] = "readable"
    state["count"] = len(data)
    return data, state


class InputOperator:
    def __init__(self, config):
        self.resources = config.get("resources", {})
        self.report = {"mode": None, "files": {}, "resources": {}, "modules": {}}
        self.skipped = {}

    # --- 1. входы и справочники ---
    def load_inputs(self, phenotype_path, variants_path, vcf_path):
        phenotypes, st = _try_load(phenotype_path, load_phenotypes)
        self.report["files"]["phenotypes"] = st

        if vcf_path:
            variants, st = _try_load(vcf_path, load_vcf)
            self.report["files"]["vcf"] = st
        else:
            variants, st = _try_load(variants_path, load_variant_table)
            self.report["files"]["variants_table"] = st

        for key in RESOURCE_KEYS:
            self.report["resources"][key] = _file_state(self.resources.get(key))

        if phenotypes and variants:
            mode = "full"
        elif phenotypes:
            mode = "phenotype_only"
        elif variants:
            mode = "variants_only"
        else:
            mode = "no_input"
        self.report["mode"] = mode
        return phenotypes, variants, mode

    def _res_ok(self, key):
        return self.report["resources"][key]["status"] == "found"

    # --- 2. какие модули можно запустить ---
    def build_analyzers(self, mode):
        analyzers, self.skipped = [], {}
        has_pheno = mode in ("full", "phenotype_only")
        has_vars = mode in ("full", "variants_only")
        g2p = self.resources.get("genes_to_phenotype_txt")

        if not has_pheno:
            self.skipped["phenotype"] = "нет фенотипов"
        elif not self._res_ok("genes_to_phenotype_txt"):
            self.skipped["phenotype"] = "нет справочника genes_to_phenotype"
        elif mode == "phenotype_only":
            analyzers.append(PhenotypeOnlyAnalyzer(g2p))
        else:
            analyzers.append(HPOAnalyzer(g2p))

        if not has_vars:
            self.skipped["pathogenicity"] = "нет вариантов"
            self.skipped["inheritance"] = "нет вариантов"
        else:
            if self._res_ok("clinvar_vcf"):
                analyzers.append(PathogenicityAnalyzer(self.resources.get("clinvar_vcf")))
            else:
                self.skipped["pathogenicity"] = "нет базы ClinVar"
            analyzers.append(InheritanceAnalyzer())
        return analyzers

    # --- 3. статусы модулей после запуска ---
    def module_statuses(self, analyzers, errors):
        status = {n: {"status": "skipped", "reason": r} for n, r in self.skipped.items()}
        for a in analyzers:
            if a.name in errors:
                status[a.name] = {"status": "failed", "error": errors[a.name]}
            else:
                status[a.name] = {"status": "ok"}
        self.report["modules"] = status
        return status

    def warnings(self):
        out = []
        for key, st in self.report["files"].items():
            if st["status"] == "missing":
                out.append(f"Файл не найден: {key} ({st['path']})")
            elif st["status"] == "unreadable":
                out.append(f"Файл не прочитан: {key}: {st['error']}")
        for key, st in self.report["resources"].items():
            if st["status"] != "found":
                out.append(f"Справочник недоступен: {key}")
        for name, st in self.report["modules"].items():
            if st["status"] == "skipped":
                out.append(f"Модуль {name} пропущен: {st['reason']}")
        return out

    def format_report(self):
        lines = ["=== Оператор ===", f"Режим: {self.report['mode']}"]
        for section, title in (("files", "Входные файлы"), ("resources", "Справочники"), ("modules", "Модули")):
            if not self.report[section]:
                continue
            lines.append(f"{title}:")
            for key, st in self.report[section].items():
                extra = st.get("error") or st.get("reason") or st.get("path") or ""
                lines.append(f"  {key}: {STATUS_RU.get(st['status'], st['status'])} {extra}".rstrip())
        return "\n".join(lines)