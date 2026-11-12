import os
from core.exceptions import ModelExecutionError
from .prompt_builder import build_factual_prompt

def explain(result, config):
    fallback = _fallback_text(result)
    settings = config.get("analysis", {})
    if not settings.get("use_openai", False):
        return fallback
    key = os.getenv(settings.get("openai_api_key_env", "OPENAI_API_KEY"))
    if not key:
        return fallback + "\n\nLLM не запущена: переменная API-ключа не задана."
    try:
        from openai import OpenAI
        client = OpenAI(api_key=key)
        prompt = build_factual_prompt(result.model_dump(mode="json"))
        response = client.chat.completions.create(
            model=settings.get("openai_model", "gpt-4o-mini"),
            messages=[{"role":"system","content":"Ты аккуратный научный ассистент по интерпретации геномных данных."},
                      {"role":"user","content":prompt}],
            temperature=0.1
        )
        return response.choices[0].message.content or fallback
    except Exception as e:
        return fallback + f"\n\nLLM fallback: запрос не выполнен ({type(e).__name__})."

def _fallback_text(result):
    if not result.candidates:
        return ("Кандидаты не ранжированы: в данных нет генов, связанных с вариантами. "
                "Это не исключает генетическую причину; проверьте аннотацию и полноту входных данных.")
    lines = ["Автоматический исследовательский отчёт (не клинический диагноз).",
             "Рейтинг основан только на доступных входных аннотациях и включённых локальных ресурсах."]
    for i, c in enumerate(result.candidates[:10], 1):
        lines.append(f"{i}. {c.gene}: интегральный балл {c.score:.3f}.")
        lines.extend(f"   - {e.statement}" for e in c.evidence)
        lines.extend(f"   Ограничение: {x}" for x in c.caveats)
    lines.append("Отсутствующая аннотация трактуется как неизвестность, а не как доброкачественность.")
    return "\n".join(lines)
