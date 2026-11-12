from pathlib import Path
from typing import Any
import yaml
from .exceptions import InvalidInputError


def load_config(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise InvalidInputError(f"Config file not found: {p}")
    with p.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise InvalidInputError("config.yaml must contain a YAML mapping")
    return data
