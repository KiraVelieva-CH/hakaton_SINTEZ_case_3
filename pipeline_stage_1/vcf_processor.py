from core.models import Variant

def qc_variants(variants: list[Variant], config: dict) -> tuple[list[Variant], list[str]]:
    t = config.get("thresholds", {})
    kept, warnings = [], []
    for v in variants:
        reasons = []
        if v.qual is not None and v.qual < float(t.get("min_variant_quality", 0)): reasons.append("low_QUAL")
        if v.depth is not None and v.depth < int(t.get("min_depth", 0)): reasons.append("low_DP")
        if v.population_af is not None and v.population_af > float(t.get("max_population_af", 1)): reasons.append("common_variant")
        if v.allele_balance is not None and not (float(t.get("min_allele_balance", 0)) <= v.allele_balance <= float(t.get("max_allele_balance", 1))):
            reasons.append("allele_balance_out_of_range")
        v.filters.extend(reasons)
        if reasons:
            warnings.append(f"{v.variant_id}: QC flags: {', '.join(reasons)}")
        # Retain flagged variants for auditability; ranking penalizes them rather than silently deleting.
        kept.append(v)
    return kept, warnings
