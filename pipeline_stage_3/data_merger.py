def merge_block_results(block_results):
    genes = {}
    for block, payload in block_results.items():
        for gene, data in payload.get("genes", {}).items():
            genes.setdefault(gene, {})[block] = data
        if block == "pathogenicity":
            for vid, data in payload.get("variants", {}).items():
                genes.setdefault("_variants", {})[vid] = data
    return genes
