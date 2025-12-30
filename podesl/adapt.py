from __future__ import annotations

import copy
from typing import Any, Callable, Dict

__all__ = ["adapt_loop"]


def adapt_loop(
    initial_model: Dict[str, Any],
    solver_func: Callable[[Dict[str, Any]], Dict[str, Any]],
    error_estimator: Callable[[Dict[str, Any], Dict[str, Any], Callable[[Dict[str, Any]], float]], float],
    refinment_strategy: Callable[[Dict[str, Any], Dict[str, Any], float], Dict[str, Any]],
    target_functional: Callable[[Dict[str, Any]], float],
    max_refine: int = 3,
    tol: float = 1e-3,
) -> Dict[str, Any]:
    model = copy.deepcopy(initial_model)
    history = {"error": [], "functional": [], "models": []}
    result = None
    for _ in range(max_refine):
        result = solver_func(model)
        J = target_functional(result)
        err = error_estimator(model, result, target_functional)
        history["functional"].append(J)
        history["error"].append(err)
        history["models"].append(model)
        if err < tol:
            break
        model = refinment_strategy(model, result, err)
    return {"final_model": model, "final_result": result, "history": history}
