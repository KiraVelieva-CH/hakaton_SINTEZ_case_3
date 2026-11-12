from __future__ import annotations
from core.models import Phenotype

def extract_phenotypes(text: str, llm_client=None) -> list[Phenotype]:
    """Conservative baseline: split supplied clinical text; optional LLM hook is intentionally not implicit."""
    if llm_client:
        return llm_client.extract_phenotypes(text)
    chunks = [x.strip(" \t-•.;") for x in text.replace("\r", "\n").split("\n") if x.strip()]
    if len(chunks) == 1:
        # Keep the original sentence intact; do not hallucinate individual HPO terms.
        return [Phenotype(term=chunks[0], source_text=chunks[0])]
    return [Phenotype(term=x, source_text=x) for x in chunks]
