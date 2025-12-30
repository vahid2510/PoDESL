from __future__ import annotations

from typing import Any, Callable, Dict, List

import numpy as np

__all__ = ["calibrate_parameters"]


def _evaluate_cost(
    params: np.ndarray,
    param_names: List[str],
    model_builder: Callable[[Dict[str, Any]], Dict[str, Any]],
    solver_func: Callable[[Dict[str, Any]], Dict[str, Any]],
    data_x: np.ndarray,
    data_y: np.ndarray,
    obs_extractor: Callable[[Dict[str, Any], np.ndarray], np.ndarray],
) -> float:
    param_dict = dict(zip(param_names, params))
    model = model_builder(param_dict)
    result = solver_func(model)
    pred = np.asarray(obs_extractor(result, data_x), dtype=float)
    residual = data_y - pred
    return 0.5 * float(np.dot(residual, residual))


def _finite_diff_cost_grad(
    params: np.ndarray,
    *args,
    eps: float = 1e-6,
) -> np.ndarray:
    grad = np.zeros_like(params)
    for i in range(params.size):
        step = eps * max(1.0, abs(params[i]))
        params_f = params.copy()
        params_b = params.copy()
        params_f[i] += step
        params_b[i] -= step
        cost_f = _evaluate_cost(params_f, *args)
        cost_b = _evaluate_cost(params_b, *args)
        grad[i] = (cost_f - cost_b) / (2.0 * step)
    return grad


def calibrate_parameters(
    param_names: List[str],
    param_init: np.ndarray,
    model_builder: Callable[[Dict[str, Any]], Dict[str, Any]],
    solver_func: Callable[[Dict[str, Any]], Dict[str, Any]],
    data_x: np.ndarray,
    data_y: np.ndarray,
    obs_extractor: Callable[[Dict[str, Any], np.ndarray], np.ndarray],
    max_iter: int = 50,
    tol: float = 1e-6,
    step_size: float = 1e-2,
) -> Dict[str, Any]:
    params = np.asarray(param_init, dtype=float).copy()
    history = {"params": [], "cost": []}
    converged = False
    args = (param_names, model_builder, solver_func, data_x, data_y, obs_extractor)

    for iteration in range(1, max_iter + 1):
        cost = _evaluate_cost(params, *args)
        history["params"].append(params.copy())
        history["cost"].append(cost)
        grad = _finite_diff_cost_grad(params, *args)
        grad_norm = np.linalg.norm(grad)
        if grad_norm < tol:
            converged = True
            break
        params -= step_size * grad

    return {
        "param_names": param_names,
        "param_opt": params,
        "history": history,
        "converged": converged,
        "iterations": iteration,
    }
