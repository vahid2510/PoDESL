from __future__ import annotations

import copy
from typing import Any, Callable, Dict, Iterable, List

import numpy as np

from ..dispatcher import dispatch
from ..uncertainty import mc_propagate
from ..optimize import optimize_design
from ..scenario import apply_modifications, run_scenarios

__all__ = [
    "solve_study_montecarlo",
    "solve_study_optimize",
    "solve_study_scenarios",
]


def _require(env: Dict[str, Any], key: str):
    if key not in env:
        raise KeyError(f"Input '{key}' required for study driver")
    return env[key]


def _solver_from_child(child_problem: str) -> Callable[[Dict[str, Any]], Dict[str, Any]]:
    def _solver(model_env: Dict[str, Any]) -> Dict[str, Any]:
        return dispatch({"title": child_problem}, model_env)

    return _solver


def solve_study_montecarlo(env: Dict[str, Any]) -> Dict[str, Any]:
    child_problem = _require(env, "child_problem")
    base_env = copy.deepcopy(_require(env, "base_env"))
    param_defs = copy.deepcopy(_require(env, "param_defs"))
    outputs = list(env.get("outputs") or [])
    nsamples = int(env.get("nsamples", 1))
    rng_seed = env.get("rng_seed")

    solver = _solver_from_child(child_problem)
    stats = mc_propagate(
        solver_func=solver,
        base_model=base_env,
        param_defs=param_defs,
        outputs=outputs,
        nsamples=nsamples,
        rng_seed=int(rng_seed) if rng_seed is not None else None,
    )
    return stats


def _design_env(base_env: Dict[str, Any], paths: Iterable[str], values: Iterable[float]) -> Dict[str, Any]:
    mods = []
    for path, value in zip(paths, values):
        mods.append({"path": path, "set": float(value)})
    return apply_modifications(base_env, mods)


def solve_study_optimize(env: Dict[str, Any]) -> Dict[str, Any]:
    child_problem = _require(env, "child_problem")
    base_env = copy.deepcopy(_require(env, "base_env"))
    design_specs = env.get("design_vars") or []
    if not design_specs:
        raise ValueError("design_vars must be provided")
    design_paths = [spec["path"] if isinstance(spec, dict) else str(spec) for spec in design_specs]
    x0 = np.asarray(env.get("x0", [0.0] * len(design_paths)), dtype=float)
    if x0.size != len(design_paths):
        raise ValueError("Length of x0 must match design_vars")
    expr = env.get("objective_expr")
    if not expr:
        raise ValueError("objective_expr is required")

    solver = _solver_from_child(child_problem)

    def f_obj(x: np.ndarray) -> float:
        model_env = _design_env(base_env, design_paths, x)
        result = solver(model_env)
        local_ctx = {**result, "np": np}
        value = eval(expr, {"__builtins__": {}}, local_ctx)
        return float(value)

    optim = optimize_design(
        f_obj,
        x0,
        max_iter=int(env.get("max_iter", 100)),
        step_size=float(env.get("step_size", 1e-2)),
        tol=float(env.get("tol", 1e-6)),
        verbose=bool(env.get("verbose", False)),
    )
    model_opt = _design_env(base_env, design_paths, optim["x_opt"])
    result_opt = solver(model_opt)
    return {
        "x_opt": optim["x_opt"],
        "f_opt": optim["f_opt"],
        "history": optim["history"],
        "converged": optim["converged"],
        "iterations": optim["iterations"],
        "model": model_opt,
        "result": result_opt,
    }


def _metric_functions(specs: Iterable[Any]) -> Dict[str, Callable[[Dict[str, Any]], Any]]:
    metrics: Dict[str, Callable[[Dict[str, Any]], Any]] = {}
    for spec in specs or []:
        if isinstance(spec, str):
            name = spec
            expr = spec
        else:
            name = spec.get("name")
            expr = spec.get("expr", name)
        if not name or not expr:
            continue

        def _make(expr_str: str) -> Callable[[Dict[str, Any]], Any]:
            def _extract(res: Dict[str, Any]) -> Any:
                ctx = {**res, "np": np}
                return eval(expr_str, {"__builtins__": {}}, ctx)

            return _extract

        metrics[name] = _make(expr)
    return metrics


def solve_study_scenarios(env: Dict[str, Any]) -> Dict[str, Any]:
    child_problem = _require(env, "child_problem")
    base_env = copy.deepcopy(_require(env, "base_env"))
    scenarios = copy.deepcopy(env.get("scenarios") or [])
    if not scenarios:
        raise ValueError("scenarios list required")
    metric_specs = env.get("metrics") or []
    metrics = _metric_functions(metric_specs)
    solver = _solver_from_child(child_problem)
    summary = run_scenarios(base_env, scenarios, solver, metrics)
    return summary
