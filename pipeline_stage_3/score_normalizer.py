def minmax(values):
    present = [v for v in values if v is not None]
    if not present: return [None for _ in values]
    lo, hi = min(present), max(present)
    if hi == lo: return [0.5 if v is not None else None for v in values]
    return [((v-lo)/(hi-lo) if v is not None else None) for v in values]
