from __future__ import annotations

from typing import Dict

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu


def _uniform_mesh(L: float, nel: int) -> np.ndarray:
    if nel <= 0:
        raise ValueError("nel must be positive.")
    return np.linspace(0.0, L, nel + 1, dtype=float)


def _assemble_conduction(L: float, A: float, k: float, nel: int):
    x = _uniform_mesh(L, nel)
    h = L / nel
    K = lil_matrix((nel + 1, nel + 1), dtype=float)
    ke = (k * A / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
    for e in range(nel):
        dofs = (e, e + 1)
        for i in range(2):
            for j in range(2):
                K[dofs[i], dofs[j]] += ke[i, j]
    return x, K


def _apply_dirichlet(K: lil_matrix, F: np.ndarray, node: int, value: float):
    F -= K[:, node].toarray().ravel() * value
    K.rows[node] = [node]
    K.data[node] = [1.0]
    for r in range(len(F)):
        if r == node:
            continue
        row = K.rows[r]
        data = K.data[r]
        if node in row:
            idx = row.index(node)
            row.pop(idx)
            data.pop(idx)
    F[node] = value


def solve_thermal_1d_steady(params: Dict[str, float | int]) -> dict:
    L = float(params["L"])
    A = float(params.get("A", 1.0))
    k = float(params["k"])
    nel = int(params.get("nel", 10))
    q_gen = float(params.get("q_gen", 0.0))

    x, K = _assemble_conduction(L, A, k, nel)
    F = np.zeros(nel + 1, dtype=float)
    F += q_gen * A * (L / nel) / 2.0

    if params.get("T_left") is not None:
        _apply_dirichlet(K, F, 0, float(params["T_left"]))
    if params.get("T_right") is not None:
        _apply_dirichlet(K, F, nel, float(params["T_right"]))

    T = splu(csc_matrix(K)).solve(F)
    return {"x": x, "T": T}


def solve_thermal_1d_transient(params: Dict[str, float | int]) -> dict:
    L = float(params["L"])
    A = float(params.get("A", 1.0))
    k = float(params["k"])
    rho = float(params["rho"])
    c = float(params["c"])
    nel = int(params.get("nel", 10))
    dt = float(params["dt"])
    t_end = float(params["t_end"])
    q_gen = float(params.get("q_gen", 0.0))

    x, K = _assemble_conduction(L, A, k, nel)
    h = L / nel
    C = lil_matrix((nel + 1, nel + 1), dtype=float)
    ce = (rho * c * A * h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
    for e in range(nel):
        dofs = (e, e + 1)
        for i in range(2):
            for j in range(2):
                C[dofs[i], dofs[j]] += ce[i, j]

    T0 = params.get("T0")
    if T0 is None:
        Tn = np.full(nel + 1, float(params.get("T_left", params.get("T_right", 300.0))))
    else:
        arr = np.asarray(T0, dtype=float)
        if arr.size == 1:
            Tn = np.full(nel + 1, float(arr))
        else:
            Tn = arr.reshape(nel + 1)

    A_sys = (C / dt + K).tolil()
    B = C / dt

    timeline = np.linspace(0.0, t_end, int(round(t_end / dt)) + 1)
    Thist = np.zeros((timeline.size, nel + 1), dtype=float)
    Thist[0, :] = Tn

    rhs_base = np.zeros(nel + 1, dtype=float)
    rhs_base += q_gen * A * h / 2.0

    for step in range(1, timeline.size):
        rhs = B @ Tn + rhs_base
        system = A_sys.copy()
        if params.get("T_left") is not None:
            _apply_dirichlet(system, rhs, 0, float(params["T_left"]))
        if params.get("T_right") is not None:
            _apply_dirichlet(system, rhs, nel, float(params["T_right"]))
        Tn1 = splu(csc_matrix(system)).solve(rhs)
        Thist[step, :] = Tn1
        Tn = Tn1

    return {"x": x, "time": timeline, "T": Thist}
