from __future__ import annotations

import copy
from typing import Any, Callable, Dict, List

import numpy as np

__all__ = ["apply_modifications", "run_scenarios"]


def _parse_path(path: str) -> List[Any]:
    segments: List[Any] = []
    token = ""
    i = 0
    while i < len(path):
        ch = path[i]
        if ch == '.':
            if token:
                segments.append(token)
                token = ""
            i += 1
        elif ch == '[':
            if token:
                segments.append(token)
                token = ""
            j = path.index(']', i)
            idx = int(path[i + 1 : j])
            segments.append(idx)
            i = j + 1
        else:
            token += ch
            i += 1
    if token:
        segments.append(token)
    return segments


def apply_modifications(model: Dict[str, Any], modifications: List[Dict[str, Any]]) -> Dict[str, Any]:
    new_model = copy.deepcopy(model)
    for mod in modifications or []:
        path = mod.get("path")
        if not path:
            continue
        segments = _parse_path(path)
        if not segments:
            continue
        parent = new_model
        for seg in segments[:-1]:
            parent = parent[seg]
        key = segments[-1]
        if "set" in mod:
            parent[key] = mod["set"]
        elif "scale" in mod:
            parent[key] = parent[key] * mod["scale"]
        elif "add" in mod:
            parent[key] = parent[key] + mod["add"]
    return new_model


def run_scenarios(
    base_model: Dict[str, Any],
    scenarios: List[Dict[str, Any]],
    solver_func: Callable[[Dict[str, Any]], Dict[str, Any]],
    metrics: Dict[str, Callable[[Dict[str, Any]], Any]],
) -> Dict[str, Any]:
    outputs = []
    for scenario in scenarios or []:
        model = apply_modifications(base_model, scenario.get("modifications", []))
        entry = {"name": scenario.get("name", "scenario"), "metrics": {}, "result": None}
        try:
            result = solver_func(model)
        except Exception:
            metric_values = {name: np.nan for name in metrics}
        else:
            entry["result"] = result
            metric_values = {}
            for name, extractor in (metrics or {}).items():
                try:
                    metric_values[name] = extractor(result)
                except Exception:
                    metric_values[name] = np.nan
        entry["metrics"] = metric_values
        outputs.append(entry)
    return {"scenarios": outputs}
