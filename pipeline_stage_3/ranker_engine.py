from core.models import Candidate, EvidenceItem

MODE_RU = {
    "AD": "AD (аутосомно-доминантный: достаточно одной изменённой копии гена)",
    "AR": "AR (аутосомно-рецессивный: нужны обе копии гена изменены)",
    "UNKNOWN": "неизвестен (гена нет в справочнике режимов, нейтральный скор 0.5)",
}
ZYG_RU = {
    "Гетерозигота": "гетерозигота (изменена одна копия)",
    "Гомозигота": "гомозигота (изменены обе копии)",
    "Референс": "референс (изменённых копий нет)",
}
COMP_RU = {"phenotype": "фенотип", "pathogenicity": "патогенность",
           "inheritance": "наследование", "exomiser": "Exomiser"}


def rank_candidates(context, block_results, config):
    weights = config.get("weights", {})
    wt = lambda k: float(weights.get(k, 0))
    n_present = sum(1 for p in context["phenotypes"] if p.present)

    variants_by_gene = {}
    for v in context["variants"]:
        if v.gene:
            variants_by_gene.setdefault(v.gene, []).append(v)

    genes = set(variants_by_gene)
    for block in ("phenotype", "inheritance"):
        genes |= set(block_results.get(block, {}).get("genes", {}))

    candidates = []
    for gene in genes:
        ph_data = block_results.get("phenotype", {}).get("genes", {}).get(gene, {})
        ph = ph_data.get("score")
        ph_matched = ph_data.get("matched", 0)
        # всего симптомов пациента: восстанавливаем из оценки (matched / score), иначе число present-признаков
        ph_total = round(ph_matched / ph) if (ph and ph_matched) else n_present

        inh_data = block_results.get("inheritance", {}).get("genes", {}).get(gene, {})
        inh = inh_data.get("score")
        zygosity = inh_data.get("zygosity", "неизвестно")
        mode = inh_data.get("mode", "UNKNOWN")

        exo = None
        gene_vars = variants_by_gene.get(gene, [])

        best_variant_id = gene_vars[0].variant_id if gene_vars else None
        path_scores, best_p = [], None
        for v in gene_vars:
            p_score = block_results.get("pathogenicity", {}).get("variants", {}).get(v.variant_id, {}).get("score")
            if p_score is not None:
                path_scores.append(p_score)
                if best_p is None or p_score > best_p:
                    best_p, best_variant_id = p_score, v.variant_id
        path = max(path_scores) if path_scores else None

        components = {"phenotype": ph, "pathogenicity": path, "inheritance": inh, "exomiser": exo}
        used, missing = [], []
        weighted, total_w = 0.0, 0.0
        for k, val in components.items():
            w = wt(k)
            if w <= 0:
                continue
            if val is None:
                missing.append(k)
            else:
                v_clamped = max(0.0, min(1.0, float(val)))
                weighted += v_clamped * w
                total_w += w
                used.append((k, v_clamped, w))

        if total_w:
            score, score_status = weighted / total_w, "calculated"
        else:
            score, score_status = 0.0, "not_calculated"

        evidence, caveats = [], []
        if ph is not None:
            evidence.append(EvidenceItem(
                source="Фенотип",
                statement=f"Совпадение по HPO: {ph:.3f} (найдено {ph_matched} из {ph_total} симптомов пациента), вес {wt('phenotype'):.2f}"))
        if inh is not None:
            evidence.append(EvidenceItem(
                source="Наследование",
                statement=f"Режим гена: {MODE_RU.get(mode, mode)}. Генотип пациента: {ZYG_RU.get(zygosity, zygosity)}. Скор: {inh:.3f}, вес {wt('inheritance'):.2f}"))
        if path is not None:
            evidence.append(EvidenceItem(
                source="Патогенность",
                statement=f"Функциональный скор (ClinVar): {path:.3f}, вес {wt('pathogenicity'):.2f}"))
        if used:
            parts = " + ".join(f"{val:.3f}×{w:.2f}" for _, val, w in used)
            text = f"Расчёт балла: ({parts}) / {total_w:.2f} = {score:.3f}"
            if missing:
                text += "; без данных, исключено из расчёта: " + ", ".join(COMP_RU[k] for k in missing)
            evidence.append(EvidenceItem(source="Расчёт", statement=text))

        if score_status == "not_calculated":
            caveats.append("Итоговый балл не рассчитан: нет ни одной оценки.")
        if not gene_vars:
            caveats.append("Нет вариантов: конкретный вариант не предсказан, только список генов.")
        if any(v.filters for v in gene_vars):
            caveats.append("Есть варианты с флагами QC (низкое качество).")
        if ph is None:
            caveats.append("Нет фенотипической ассоциации.")
        if path is None and gene_vars:
            caveats.append("Нет аннотации патогенности (NA в таблице).")

        candidates.append(Candidate(
            gene=gene, variant_id=best_variant_id, score=score, score_status=score_status,
            component_scores=components, evidence=evidence, caveats=caveats))

    candidates.sort(key=lambda c: c.score, reverse=True)
    return candidates