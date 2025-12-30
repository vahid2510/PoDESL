from __future__ import annotations

import numpy as np

from .bc_utils import PLANAR_UV_DOF_MAP, dirichlet_vector
from ..linalg_utils import safe_solve


def _elem_k(EA, L, c, s):
    k = np.array(
        [
            [c * c, c * s, -c * c, -c * s],
            [c * s, s * s, -c * s, -s * s],
            [-c * c, -c * s, c * c, c * s],
            [-c * s, -s * s, c * s, s * s],
        ],
        float,
    )
    return (EA / L) * k


def _property_array(env, key, count):
    values = env.get(key)
    if values is None:
        raise KeyError(f"Input '{key}' required for truss analysis.")
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        return np.full(count, float(arr))
    arr = arr.reshape(-1)
    if arr.size == 1:
        return np.full(count, float(arr[0]))
    if arr.size != count:
        raise ValueError(f"Property '{key}' must have length {count}, got {arr.size}.")
    return arr


def solve_truss2d_static(env):
    nodes = np.asarray(env.get("nodes"), float)
    elems = np.asarray(env.get("elems"), int)
    ne = elems.shape[0]
    nn = nodes.shape[0]
    ndof = 2 * nn

    E = _property_array(env, "E", ne)
    A = _property_array(env, "A", ne)

    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)
    lengths = np.zeros(ne)

    for ee, (n1, n2) in enumerate(elems):
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]
        L = float(np.hypot(x2 - x1, y2 - y1))
        if L <= 0:
            raise ValueError("Zero-length truss element")
        lengths[ee] = L
        c = (x2 - x1) / L
        s = (y2 - y1) / L
        ke = _elem_k(E[ee] * A[ee], L, c, s)
        dofs = [2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1]
        for i, I in enumerate(dofs):
            for j, J in enumerate(dofs):
                K[I, J] += ke[i, j]

    load_entries = list(env.get("loads", []) or []) + list(env.get("point_loads", []) or [])
    for ent in load_entries:
        n = int(ent[0])
        fx = float(ent[1])
        fy = float(ent[2])
        F[2 * n] += fx
        F[2 * n + 1] += fy

    fixed_idx, fixed_vals = dirichlet_vector(
        env.get("fix", []),
        ndof=ndof,
        ndof_per_node=2,
        dof_map=PLANAR_UV_DOF_MAP,
    )
    free = np.setdiff1d(np.arange(ndof), fixed_idx)

    U = np.zeros(ndof)
    if free.size > 0:
        rhs = F[free]
        if fixed_idx.size:
            rhs -= K[np.ix_(free, fixed_idx)] @ fixed_vals[fixed_idx]
        U[free] = safe_solve(K[np.ix_(free, free)], rhs)
    if fixed_idx.size:
        U[fixed_idx] = fixed_vals[fixed_idx]

    reactions = K @ U - F

    member_forces = np.zeros(ne)
    for ee, (n1, n2) in enumerate(elems):
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]
        _, _, L, T = _elem_transform(x1, y1, x2, y2)
        dofs = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1])
        u_local = T @ U[dofs]
        member_forces[ee] = (E[ee] * A[ee] / L) * (u_local[1] - u_local[0])

    disp = U.reshape(nn, 2)

    return {
        "U": U,
        "disp": disp,
        "R": reactions,
        "reactions": reactions,
        "member_forces": member_forces,
        "lengths": lengths,
        "K": K,
        "nodes": nodes,
        "elems": elems,
    }


def _elem_transform(x1, y1, x2, y2):
    L = float(np.hypot(x2 - x1, y2 - y1))
    if L <= 0:
        raise ValueError("Zero-length truss element")
    c = (x2 - x1) / L
    s = (y2 - y1) / L
    T = np.array([[c, s, 0, 0], [0, 0, c, s]], float)
    return c, s, L, T
