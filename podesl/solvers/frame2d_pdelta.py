from __future__ import annotations

import numpy as np

from .bc_utils import FRAME2D_DOF_MAP, dirichlet_vector
from .input_utils import merged_dirichlet, merged_point_loads
from ..linalg_utils import safe_solve


def _T(c, s):
    R = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]], float)
    T = np.zeros((6, 6), float)
    T[:3, :3] = R
    T[3:, 3:] = R
    return T


def _geom(x1, y1, x2, y2):
    L = float(np.hypot(x2 - x1, y2 - y1))
    if L <= 0:
        raise ValueError("Zero-length frame element")
    c = (x2 - x1) / L
    s = (y2 - y1) / L
    return c, s, L


def _k_local(EA, EI, L):
    L2 = L * L
    L3 = L2 * L
    k = np.zeros((6, 6), float)
    k[0, 0] = k[3, 3] = EA / L
    k[0, 3] = k[3, 0] = -EA / L
    k[1, 1] = k[4, 4] = 12 * EI / L3
    k[1, 4] = k[4, 1] = -12 * EI / L3
    k[1, 2] = k[2, 1] = 6 * EI / L2
    k[1, 5] = k[5, 1] = 6 * EI / L2
    k[2, 2] = k[5, 5] = 4 * EI / L
    k[2, 5] = k[5, 2] = 2 * EI / L
    k[2, 4] = k[4, 2] = -6 * EI / L2
    k[4, 5] = k[5, 4] = -6 * EI / L2
    return k


def _kg_local(N, L):
    L2 = L * L
    c1 = N / (30.0 * L)
    return c1 * np.array(
        [
            [0, 0, 0, 0, 0, 0],
            [0, 36, 3 * L, 0, -36, 3 * L],
            [0, 3 * L, 4 * L2, 0, -3 * L, -1 * L2],
            [0, 0, 0, 0, 0, 0],
            [0, -36, -3 * L, 0, 36, -3 * L],
            [0, 3 * L, -1 * L2, 0, -3 * L, 4 * L2],
        ],
        float,
    )


def _element_array(value, count, name):
    if value is None:
        raise KeyError(f"Input '{name}' required for frame analysis.")
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return np.full(count, float(arr))
    arr = arr.reshape(-1)
    if arr.size == 1:
        return np.full(count, float(arr[0]))
    if arr.size != count:
        raise ValueError(f"Property '{name}' must have length {count}, got {arr.size}.")
    return arr


def solve_frame2d_pdelta(env):
    nodes = np.asarray(env.get("nodes"), float)
    elems = np.asarray(env.get("elems"), int)
    ne = elems.shape[0]
    nn = nodes.shape[0]
    ndof = 3 * nn

    E = _element_array(env.get("E"), ne, "E")
    A = _element_array(env.get("A"), ne, "A")
    I = _element_array(env.get("I"), ne, "I")

    Nin = env.get("N", 0.0)
    N = (
        np.full(ne, float(Nin))
        if np.ndim(Nin) == 0
        else np.asarray(Nin, float).reshape(ne)
    )

    K = np.zeros((ndof, ndof), float)
    Kg = np.zeros((ndof, ndof), float)
    F = np.zeros(ndof, float)
    elem_data = []

    for ee, (n1, n2) in enumerate(elems):
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]
        c, s, L = _geom(x1, y1, x2, y2)
        kL = _k_local(E[ee] * A[ee], E[ee] * I[ee], L)
        kgL = _kg_local(N[ee], L)
        T = _T(c, s)
        kG = T.T @ kL @ T
        kgG = T.T @ kgL @ T
        dofs = np.array([3 * n1, 3 * n1 + 1, 3 * n1 + 2, 3 * n2, 3 * n2 + 1, 3 * n2 + 2])
        for i, Iglob in enumerate(dofs):
            for j, Jglob in enumerate(dofs):
                K[Iglob, Jglob] += kG[i, j]
                Kg[Iglob, Jglob] += kgG[i, j]
        elem_data.append((dofs, T, kL, kgL))

    load_entries = merged_point_loads(
        env.get("loads"),
        env.get("point_loads"),
        dof_order=("ux", "uy", "rz"),
    )
    for n, comp in load_entries:
        fx, fy = comp[0], comp[1]
        mz = comp[2] if len(comp) > 2 else 0.0
        F[3 * n + 0] += fx
        F[3 * n + 1] += fy
        F[3 * n + 2] += mz

    fix_entries = merged_dirichlet(
        env.get("fix"),
        env.get("bcs"),
        dof_map=FRAME2D_DOF_MAP,
    )
    fixed_idx, fixed_vals = dirichlet_vector(
        fix_entries,
        ndof=ndof,
        ndof_per_node=3,
        dof_map=FRAME2D_DOF_MAP,
    )
    free_idx = np.setdiff1d(np.arange(ndof), fixed_idx)

    controls = env.get("CONTROLS", {}) or {}
    load_steps = int(controls.get("load_steps", env.get("load_steps", 10)))
    max_iter = int(controls.get("max_iter", env.get("max_iter", 40)))
    tol = float(controls.get("tol", env.get("tol", 1e-6)))

    U = np.zeros(ndof)
    history = []
    for step in range(1, load_steps + 1):
        lam = step / load_steps
        F_target = lam * F
        for it in range(1, max_iter + 1):
            Ktot = K + Kg
            R = F_target - Ktot @ U
            res_norm = np.linalg.norm(R[free_idx], ord=np.inf)
            history.append(res_norm)
            if res_norm < tol:
                break
            delta = np.zeros(ndof)
            delta[free_idx] = safe_solve(Ktot[np.ix_(free_idx, free_idx)], R[free_idx])
            U += delta
            if it == max_iter:
                break

    Ktot = K + Kg
    R = Ktot @ U - F
    member_forces = np.zeros((ne, 6))
    for ee, (dofs, T, kL, kgL) in enumerate(elem_data):
        u_local = T @ U[dofs]
        member_forces[ee, :] = (kL + kgL) @ u_local

    disp = U.reshape(nn, 3)

    return {
        "U": U,
        "disp": disp,
        "R": R,
        "reactions": R,
        "K": K,
        "Kg": Kg,
        "Ktot": Ktot,
        "N": N,
        "member_forces": member_forces,
        "nodes": nodes,
        "elems": elems,
        "residual_history": history,
    }
