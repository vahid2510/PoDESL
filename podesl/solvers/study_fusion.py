from __future__ import annotations

import copy
from importlib import import_module
from typing import Any, Callable, Dict, List

import numpy as np

from ..fusion import CouplingLink, FusionEngine, Subsystem

__all__ = ["solve_study_fusion"]


def _safe_eval(expr: str, context: Dict[str, Any]) -> Any:
    env = {"__builtins__": {}, "np": np}
    env.update(context or {})
    return eval(expr, env, {})


def _resolve_callable(ref: Any) -> Callable:
    if callable(ref):
        return ref
    if not isinstance(ref, str):
        raise TypeError("Solver reference must be callable or import path string")
    if ":" in ref:
        module_name, func_name = ref.split(":", 1)
    else:
        module_name, func_name = ref.rsplit(".", 1)
    module = import_module(module_name)
    return getattr(module, func_name)


def _build_subsystems(specs: List[Dict[str, Any]]) -> List[Subsystem]:
    subsystems: List[Subsystem] = []
    for spec in specs or []:
        name = spec["name"]
        solver = _resolve_callable(spec["solver"])
        model = copy.deepcopy(spec.get("model", {}))
        subsystems.append(Subsystem(name, solver, model))
    return subsystems


def _build_couplings(specs: List[Dict[str, Any]], context: Dict[str, Any]) -> List[CouplingLink]:
    couplings: List[CouplingLink] = []
    for spec in specs or []:
        map_func = spec.get("map_func") or spec.get("map")
        if map_func is None:
            raise ValueError("Coupling specification must include 'map' or 'map_func'.")
        if isinstance(map_func, str):
            map_func = _safe_eval(map_func, context)
        if not callable(map_func):
            raise TypeError("map_func must be callable")
        couplings.append(
            CouplingLink(
                src=spec["src"],
                src_field=spec["src_field"],
                dst=spec["dst"],
                dst_field=spec["dst_field"],
                map_func=map_func,
                description=spec.get("description", ""),
            )
        )
    return couplings


def _build_time_update(spec: Any, context: Dict[str, Any]) -> Callable[[int, float, Dict[str, Subsystem]], None]:
    if spec is None:
        return lambda step, t, subsystems: None
    if callable(spec):
        return spec
    if isinstance(spec, str):
        func = _safe_eval(spec, context)
        if not callable(func):
            raise TypeError("time_update expression must evaluate to a callable")
        return func
    raise TypeError("time_update must be callable or string expression")


def solve_study_fusion(env: Dict[str, Any]) -> Dict[str, Any]:
    subspecs = env.get("subsystems")
    if not subspecs:
        raise KeyError("Fusion study requires a 'subsystems' list")
    context = env.get("context", {}) or {}
    subsystems = _build_subsystems(subspecs)
    couplings = _build_couplings(env.get("couplings", []), context)
    engine = FusionEngine(
        subsystems,
        couplings,
        mode=str(env.get("mode", "steady")).lower(),
        max_iter=int(env.get("max_iter", 10)),
        tol=float(env.get("tol", 1e-3)),
        verbose=bool(env.get("verbose", False)),
    )

    if engine.mode == "steady":
        results = engine.run_steady()
        return {
            "mode": "steady",
            "results": results,
            "residual_history": engine.residual_history,
        }

    if engine.mode == "transient":
        nsteps = int(env.get("nsteps", 1))
        dt = float(env.get("dt", 1.0))
        time_update = _build_time_update(env.get("time_update"), context)
        times, history = engine.run_transient(nsteps=nsteps, time_update=time_update, dt=dt)
        return {
            "mode": "transient",
            "times": times,
            "history": history,
            "residual_history": engine.residual_history,
        }

    raise ValueError(f"Unsupported fusion mode '{engine.mode}'")
