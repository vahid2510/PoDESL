from __future__ import annotations

from typing import Callable, Dict

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from .frame2d import solve_frame2d_static


def _loads_from_pressure(nodes: np.ndarray, pressure: Callable[[float], float], time: float, boundary: list[int]) -> list[list[float]]:
    p = float(pressure(time))
    loads = []
    for node in boundary or []:
        loads.append([node, 0.0, -p, 0.0])
    return loads


def solve_fsi2d_coupled(model: Dict[str, object]) -> Dict[str, np.ndarray]:
    structure = model["structure_model"].copy()
    pressure_func = model["pressure_time_func"]
    dt = float(model["dt"])
    t_end = float(model["t_end"])
    boundary = structure.get("pressure_nodes", [])

    nsteps = int(round(t_end / dt))
    times = np.linspace(0.0, t_end, nsteps + 1)
    history = []

    base_loads = structure.get("loads", [])
    for t in times:
        struct = structure.copy()
        struct["loads"] = base_loads + _loads_from_pressure(np.asarray(struct["nodes"], float), pressure_func, t, boundary)
        res = solve_frame2d_static(struct)
        history.append(res["disp"].copy())

    return {"time": times, "disp": np.array(history)}
