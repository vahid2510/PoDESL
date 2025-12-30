from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

__all__ = ["AssertionFailure", "check_assertions", "compute_safety_score"]


class AssertionFailure(Exception):
    pass


def _max_less_than(result: Dict[str, Any], assertion: Dict[str, Any]) -> Dict[str, Any]:
    field = assertion.get("field")
    limit = assertion.get("limit")
    message = assertion.get("message", "")
    if field not in result:
        return {"assertion": assertion, "ok": False, "value": None, "details": f"Missing field '{field}'"}
    values = np.asarray(result[field])
    max_val = float(np.max(values))
    ok = max_val < float(limit)
    details = message or f"max({field}) = {max_val:.3e}, limit = {limit:.3e}"
    return {"assertion": assertion, "ok": ok, "value": max_val, "details": details}


def _value_in_range(result: Dict[str, Any], assertion: Dict[str, Any]) -> Dict[str, Any]:
    expr = assertion.get("expr")
    min_val = assertion.get("min", -np.inf)
    max_val = assertion.get("max", np.inf)
    message = assertion.get("message", "")
    context = {**result, "np": np}
    try:
        value = float(eval(expr, {"__builtins__": {}}, context))
    except Exception:
        return {"assertion": assertion, "ok": False, "value": None, "details": f"Failed to eval '{expr}'"}
    ok = (value >= min_val) and (value <= max_val)
    details = message or f"{value:.3e} not in [{min_val}, {max_val}]"
    return {"assertion": assertion, "ok": ok, "value": value, "details": details}


_HANDLERS = {
    "max_less_than": _max_less_than,
    "value_in_range": _value_in_range,
}


def check_assertions(result: Dict[str, Any], assertions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    reports: List[Dict[str, Any]] = []
    for assertion in assertions or []:
        handler = _HANDLERS.get(assertion.get("type"))
        if handler is None:
            reports.append({"assertion": assertion, "ok": False, "value": None, "details": "Unsupported assertion type"})
            continue
        reports.append(handler(result, assertion))
    return reports


def compute_safety_score(assertion_results: List[Dict[str, Any]]) -> float:
    score = 1.0
    for item in assertion_results or []:
        if not item.get("ok", False):
            score -= 0.3
    return float(np.clip(score, 0.0, 1.0))
