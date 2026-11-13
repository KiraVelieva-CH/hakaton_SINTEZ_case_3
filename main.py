#main.py
from __future__ import annotations
import argparse, logging
from pathlib import Path
from core.config_loader import load_config
from core.models import AnalysisResult
from io_handlers.data_loader import load_phenotypes, load_variant_table, load_vcf
from io_handlers.data_export import export_result
from pipeline_stage_0.input_operator import InputOperator
from pipeline_stage_1.hpo_standardizer import HPOStandardizer
from pipeline_stage_1.vcf_processor import qc_variants
from pipeline_stage_1.vep_annotator import VEPAnnotator
from pipeline_stage_2.pheno_analyzer import PhenotypeAnalyzer
from pipeline_stage_2.patho_analyzer import PathogenicityAnalyzer
from pipeline_stage_2.inherit_analyzer import InheritanceAnalyzer
from pipeline_stage_2.hpo_analyzer import HPOAnalyzer
from pipeline_stage_2.parallel_runner import run_parallel
from pipeline_stage_3.ranker_engine import rank_candidates
from pipeline_stage_3.report_generator import build_gene_ranking
from pipeline_stage_4.consistency_checker import check_consistency
from pipeline_stage_4.llm_explainer import explain


def run_pipeline(config, phenotype_path=None, variants_path=None, vcf_path=None, case_id=None):
    inp, resources = config.get("input", {}), config.get("resources", {})
    phenotype_path = phenotype_path or inp.get("phenotype_file")
    variants_path = variants_path or inp.get("variants_file")
    vcf_path = vcf_path or inp.get("vcf_file")
    default_id = Path(phenotype_path).stem.replace("_phenotype", "") if phenotype_path else "case"
    case_id = case_id or inp.get("case_id") or default_id
    # Оператор: проверка и чтение входов
    operator = InputOperator(config)
    phenotypes, variants, mode = operator.load_inputs(phenotype_path, variants_path, vcf_path)
    if mode == "no_input":
        print(operator.format_report())
        raise ValueError("Нет ни читаемых фенотипов, ни вариантов: анализировать нечего")

    warnings = []
    if phenotypes:
        phenotypes = HPOStandardizer(resources.get("hpo_terms_tsv")).standardize(phenotypes)
    if variants:
        variants = VEPAnnotator(config).annotate(variants)
        n_no_gene = sum(1 for v in variants if not v.gene)
        print(f"Аннотация VEP: ген определён у {len(variants) - n_no_gene} из {len(variants)} вариантов")
        if n_no_gene:
            warnings.append(f"VEP: у {n_no_gene} из {len(variants)} вариантов нет гена "
                            "(аннотация недоступна или вариант вне генов), они не участвуют в ранжировании")
        operator.report["annotation"] = {"variants_total": len(variants), "without_gene": n_no_gene}
        variants, qc_warnings = qc_variants(variants, config)
        warnings.extend(qc_warnings)

    context = {"case_id": case_id, "phenotypes": phenotypes, "variants": variants,
               "genome_build": config.get("app", {}).get("genome_build", "GRCh38")}

    analyzers = operator.build_analyzers(mode)
    blocks, errors = run_parallel(analyzers, context, config.get("analysis", {}).get("workers", 4))
    warnings.extend(f"{k} analyzer failed: {v}" for k, v in errors.items())
    modules = operator.module_statuses(analyzers, errors)
    warnings.extend(operator.warnings())
    print(operator.format_report())

    candidates = rank_candidates(context, blocks, config)
    candidates = candidates[:int(config.get("analysis", {}).get("top_n", 50))]
    all_ok = all(m["status"] == "ok" for m in modules.values())
    result = AnalysisResult(
        case_id=case_id, genome_build=context["genome_build"],
        phenotypes=phenotypes, variants=variants, candidates=candidates,
        gene_ranking=build_gene_ranking(candidates, variants), block_results=blocks,
        consistency=check_consistency(blocks), warnings=warnings,
        status="completed" if all_ok else "partial",
        metadata={"analysis_errors": errors, "operator": operator.report,
                  "scoring_note": "Weighted mean over available components; missing components are excluded from denominator."})
    result.explanation = explain(result, config)
    return result


def main():
    parser = argparse.ArgumentParser(description="Genomic Detective research prototype")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--data-dir", help="Папка с кейсами (development или hidden)")
    parser.add_argument("--case-id", help="ID кейса (например, case_01)")
    parser.add_argument("--phenotypes", help="Прямой путь к JSON/TSV фенотипа")
    parser.add_argument("--vcf", help="Прямой путь к VCF")
    parser.add_argument("--output")
    args = parser.parse_args()
    
    config = load_config(args.config)
    logging.basicConfig(level=getattr(logging, config.get("app",{}).get("log_level","INFO").upper(), logging.INFO))
    
   
    pheno_path = args.phenotypes
    vcf_path = args.vcf
    case_id = args.case_id
    
    if args.data_dir and args.case_id:
        base_dir = Path(args.data_dir)
        pheno_path = base_dir / f"{args.case_id}_phenotype.json"
        vcf_path = base_dir / f"{args.case_id}.vcf"
        
        if not pheno_path.exists():
            raise FileNotFoundError(f"Не найден фенотип: {pheno_path}")
        if not vcf_path.exists():
            raise FileNotFoundError(f"Не найден VCF: {vcf_path}")

    result = run_pipeline(
        config, 
        phenotype_path=str(pheno_path) if pheno_path else None, 
        vcf_path=str(vcf_path) if vcf_path else None, 
        case_id=case_id
    )
    
    paths = export_result(result, args.output or config.get("app",{}).get("output_dir","output"), min_score_threshold=0.3)
    
    print(f"Analysis complete: {paths['json']} and {paths['predictions_tsv']}")
    print(result.explanation)

if __name__ == "__main__":
    main()
