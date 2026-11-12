def check_consistency(block_results):
    issues = []
    ph = block_results.get("phenotype", {}).get("genes", {})
    inh = block_results.get("inheritance", {}).get("genes", {})
    path = block_results.get("pathogenicity", {}).get("variants", {})
    for gene, pdata in ph.items():
        if pdata.get("score", 0) > 0.5 and gene not in inh:
            issues.append({"type":"evidence_disagreement", "gene":gene,
                           "message":"Phenotype associations exist, but no variant-level inheritance evidence was available."})
    for vid, pdata in path.items():
        sig = str(pdata.get("clinvar", {}).get("significance", "")).lower()
        score = pdata.get("score")
        if "benign" in sig and score is not None and score > 0.7:
            issues.append({"type":"annotation_conflict","variant_id":vid,
                           "message":"ClinVar text suggests benign while aggregate pathogenicity score is high; inspect source annotations."})
    return issues
