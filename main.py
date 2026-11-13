# main.py
from __future__ import annotations
import argparse, logging
from pathlib import Path
from core.config_loader import load_config
from core.models import AnalysisResult
from io_handlers.data_loader import load_phenotypes, load_variant_table, load_vcf
from io_handlers.data_export import export_result
from pipeline_stage_1.hpo_standardizer import HPOStandardizer
from pipeline_stage_1.vcf_processor import qc_variants
from pipeline_stage_1.vep_annotator import VEPAnnotator
from pipeline_stage_2.hpo_analyzer import HPOAnalyzer
from pipeline_stage_2.patho_analyzer import PathogenicityAnalyzer
from pipeline_stage_2.inherit_analyzer import InheritanceAnalyzer
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
    case_id = case_id or inp.get("case_id", "case")
    
    if not phenotype_path: 
        raise ValueError("Provide --phenotypes or input.phenotype_file")
        
    phenotypes = load_phenotypes(phenotype_path)
    variants = load_vcf(vcf_path) if vcf_path else load_variant_table(variants_path)
    
    phenotypes = HPOStandardizer(resources.get("hpo_terms_tsv")).standardize(phenotypes)
    
    
    variants = VEPAnnotator(config).annotate(variants)
    variants, warnings = qc_variants(variants, config)

    context = {"case_id": case_id, "phenotypes": phenotypes, "variants": variants,
               "genome_build": config.get("app", {}).get("genome_build", "GRCh38")}
               
    analyzers = [
        HPOAnalyzer(resources.get("genes_to_phenotype_txt")), 
        PathogenicityAnalyzer(resources.get("clinvar_vcf")),
        InheritanceAnalyzer(),
    ]
    
    blocks, errors = run_parallel(analyzers, context, config.get("analysis", {}).get("workers", 4))
    warnings.extend(f"{k} analyzer failed: {v}" for k,v in errors.items())
    
    candidates = rank_candidates(context, blocks, config)
    candidates = candidates[:int(config.get("analysis", {}).get("top_n", 50))]
    
    result = AnalysisResult(
        case_id=case_id, 
        genome_build=context["genome_build"],
        phenotypes=phenotypes, 
        variants=variants, 
        candidates=candidates,
        gene_ranking=build_gene_ranking(candidates, variants), 
        block_results=blocks,
        consistency=check_consistency(blocks), 
        warnings=warnings,
        metadata={
            "analysis_errors": errors, 
            "scoring_note": "Weighted mean over available components; missing components are excluded from denominator."
        }
    )
    
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
