import os
from core.exceptions import ModelExecutionError
from .prompt_builder import build_factual_prompt

def explain(result, config):
    
    return _fallback_text(result)

def _fallback_text(result):
    if not result.candidates:
        return "Кандидаты не найдены. Рекомендуется ответ no_causal_variant_in_vcf."
        
    lines = [
        "==========================================================",
        "ОТЧЕТ О РАНЖИРОВАНИИ ГЕНЕТИЧЕСКИХ ВАРИАНТОВ",
        "==========================================================",
        f"Кейс: {result.case_id}",
        "----------------------------------------------------------"
    ]
    
    for i, c in enumerate(result.candidates[:3], 1): 
        lines.append(f"КАНДИДАТ #{i}: Ген {c.gene} (Интегральный балл: {c.score:.3f})")
        lines.append(f"  Вариант: {c.variant_id or 'не указан'}")
        
        
        for e in c.evidence:
            lines.append(f"  [+] {e.statement}")
            
        
        for cav in c.caveats:
            lines.append(f"  [!] Ограничение: {cav}")
        lines.append("----------------------------------------------------------")
        
    lines.append("Отсутствующая аннотация трактуется как неизвестность, а не как доброкачественность.")
    return "\n".join(lines)
