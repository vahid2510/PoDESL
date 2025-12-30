from __future__ import annotations

import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve

from .bc_utils import PLANAR_UV_DOF_MAP, dirichlet_vector
from .input_utils import merged_dirichlet, merged_point_loads


def _elem_transform(x1, y1, x2, y2):
    L = float(np.hypot(x2 - x1, y2 - y1))
    if L <= 0:
        raise ValueError("Zero-length truss element")
    c = (x2 - x1) / L
    s = (y2 - y1) / L
    T = np.array([[c, s, 0, 0], [0, 0, c, s]], float)
    return c, s, L, T


def _element_property(env, keys, count):
    if isinstance(keys, str):
        keys = [keys]
    values = None
    for key in keys:
        if key in env and env[key] is not None:
            values = env[key]
            break
    if values is None:
        raise KeyError(f"Input '{keys[0]}' required for truss analysis.")
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        return np.full(count, float(arr))
    arr = arr.reshape(-1)
    if arr.size == 1:
        return np.full(count, float(arr[0]))
    if arr.size != count:
        raise ValueError(f"Property '{keys[0]}' must have length {count}, got {arr.size}.")
    return arr


def solve_truss2d(env):
    node_array = np.asarray(env.get("nodes"), float)
    elems = env.get("elems")
    if elems is None:
        elems = env.get("elements")
    elems = np.asarray(elems, int)

    nn = node_array.shape[0]
    ne = elems.shape[0]
    ndof = 2 * nn

    E = _element_property(env, ["E"], ne)
    A = _element_property(env, ["A", "areas"], ne)

    K = lil_matrix((ndof, ndof), dtype=float)
    F = np.zeros(ndof, dtype=float)
    Ls = np.zeros(ne, dtype=float)

    for e, (n1, n2) in enumerate(elems):
        x1, y1 = node_array[n1]
        x2, y2 = node_array[n2]
        c, s, L, T = _elem_transform(x1, y1, x2, y2)
        Ls[e] = L
        k_loc = (E[e] * A[e] / L) * np.array([[1, -1], [-1, 1]], float)
        k_glob = T.T @ k_loc @ T
        dofs = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1])
        for i, I in enumerate(dofs):
            for j, J in enumerate(dofs):
                K[I, J] += k_glob[i, j]

    load_entries = merged_point_loads(
        env.get("loads"),
        env.get("point_loads"),
        dof_order=("ux", "uy"),
    )
    for node, values in load_entries:
        F[2 * node] += values[0]
        F[2 * node + 1] += values[1]

    bc_entries = merged_dirichlet(
        env.get("fix"),
        env.get("bcs"),
        dof_map=PLANAR_UV_DOF_MAP,
    )
    fixed_idx, fixed_vals = dirichlet_vector(
        bc_entries,
        ndof=ndof,
        ndof_per_node=2,
        dof_map=PLANAR_UV_DOF_MAP,
    )
    free_idx = np.setdiff1d(np.arange(ndof), fixed_idx)

    U = np.zeros(ndof, dtype=float)
    K_csc = csc_matrix(K)

    if free_idx.size:
        rhs = F[free_idx].copy()
        if fixed_idx.size:
            rhs -= K_csc[free_idx][:, fixed_idx] @ fixed_vals[fixed_idx]
        U[free_idx] = spsolve(K_csc[free_idx][:, free_idx], rhs)
    if fixed_idx.size:
        U[fixed_idx] = fixed_vals[fixed_idx]

    reactions = (K @ U) - F

    axial = np.zeros(ne, dtype=float)
    for e, (n1, n2) in enumerate(elems):
        x1, y1 = node_array[n1]
        x2, y2 = node_array[n2]
        c, s, L, T = _elem_transform(x1, y1, x2, y2)
        dofs = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1])
        u_local = T @ U[dofs]
        axial[e] = (E[e] * A[e] / L) * (u_local[1] - u_local[0])

    out = {
        "u": U,
        "reactions": reactions,
        "element_forces": axial,
        "nodes": node_array,
        "elems": elems,
        "L": Ls,
    }
    out.update(
        {
            "U": U,
            "disp": U.reshape(nn, 2),
            "R": reactions,
            "member_forces": axial,
            "member_end_forces": None,
        }
    )
    return out
