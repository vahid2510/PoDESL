from __future__ import annotations

from typing import Dict

from .frame2d import solve_frame2d_static as _frame2d_static
from .frame2d_pdelta import solve_frame2d_pdelta as _frame2d_pdelta


def _build_env(model: Dict[str, object]) -> Dict[str, object]:
    env = {
        "nodes": model["nodes"],
        "elems": model["elements"],
        "E": model["E"],
        "I": model["I"],
        "A": model["A"],
        "loads": model.get("loads"),
        "point_loads": model.get("point_loads"),
        "fix": model.get("fix"),
        "bcs": model.get("bcs"),
    }
    return env


def solve_frame2d_beamcolumn_static(model: Dict[str, object]) -> Dict[str, object]:
    env = _build_env(model)
    return _frame2d_static(env)


def solve_frame2d_beamcolumn_nonlinear(model: Dict[str, object]) -> Dict[str, object]:
    env = _build_env(model)
    controls = {
        "method": "loadstep",
        "load_steps": int(model.get("load_steps", 10)),
        "max_iter": int(model.get("max_iter", 40)),
        "tol": float(model.get("tol", 1e-6)),
        "line_search": bool(model.get("line_search", True)),
        "report_history": True,
    }
    env["CONTROLS"] = controls
    env["nonlinear"] = True
    return _frame2d_pdelta(env)
