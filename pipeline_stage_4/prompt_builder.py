def build_factual_prompt(result_dict):
    # Only facts already present in the analysis are passed to the model.
    return (
        "Explain the following computational research result. Do not diagnose. "
        "Do not add genes, variants, phenotype associations, or database facts not present in the JSON. "
        "Clearly distinguish measured/input facts, heuristic scores, missing data, and uncertainty. "
        "If evidence is insufficient, say so. Return concise Russian text.\nJSON:\n"
        + __import__("json").dumps(result_dict, ensure_ascii=False, indent=2)
    )
