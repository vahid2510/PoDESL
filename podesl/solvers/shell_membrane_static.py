from __future__ import annotations

import numpy as np

from .bc_utils import PLANAR_UV_DOF_MAP, dirichlet_vector
from .input_utils import merged_dirichlet, merged_point_loads

__all__ = ["solve_shell_membrane_static"]


def _constitutive_matrix(E: float, nu: float) -> np.ndarray:
    coeff = E / (1.0 - nu * nu)
    D = coeff * np.array(
        [
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ],
        dtype=float,
    )
    return D


def _triangle_B_matrix(x1, y1, x2, y2, x3, y3) -> tuple[np.ndarray, float]:
    A2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    area = 0.5 * A2
    if area <= 0:
        raise ValueError("Degenerate or inverted triangle element in shell membrane mesh")
    beta = np.array([y2 - y3, y3 - y1, y1 - y2], dtype=float)
    gamma = np.array([x3 - x2, x1 - x3, x2 - x1], dtype=float)
    B = np.zeros((3, 6), dtype=float)
    for i in range(3):
        B[0, 2 * i] = beta[i]
        B[1, 2 * i + 1] = gamma[i]
        B[2, 2 * i] = gamma[i]
        B[2, 2 * i + 1] = beta[i]
    B /= (2.0 * area)
    return B, area


def solve_shell_membrane_static(model: dict) -> dict:
    nodes = np.asarray(model["nodes"], dtype=float)
    elements = np.asarray(model["elements"], dtype=int)
    nn = nodes.shape[0]
    ndof = 2 * nn

    t = float(model.get("thickness", 1.0))
    E = float(model["E"])
    nu = float(model.get("nu", 0.3))
    D = _constitutive_matrix(E, nu)

    K = np.zeros((ndof, ndof), dtype=float)
    F = np.zeros(ndof, dtype=float)
    elem_sigma = np.zeros((elements.shape[0], 3), dtype=float)

    for ee, (n1, n2, n3) in enumerate(elements):
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]
        x3, y3 = nodes[n3]
        B, area = _triangle_B_matrix(x1, y1, x2, y2, x3, y3)
        ke = t * area * (B.T @ D @ B)
        dofs = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n3, 2 * n3 + 1])
        for i, I in enumerate(dofs):
            for j, J in enumerate(dofs):
                K[I, J] += ke[i, j]

    load_entries = merged_point_loads(
        model.get("loads"),
        model.get("point_loads"),
        dof_order=("ux", "uy"),
    )
    for node, values in load_entries:
        F[2 * node] += values[0]
        F[2 * node + 1] += values[1]

    fix_entries = merged_dirichlet(
        model.get("fix"),
        model.get("bcs"),
        dof_map=PLANAR_UV_DOF_MAP,
    )
    fixed_idx, fixed_vals = dirichlet_vector(
        fix_entries,
        ndof=ndof,
        ndof_per_node=2,
        dof_map=PLANAR_UV_DOF_MAP,
    )
    free_idx = np.setdiff1d(np.arange(ndof), fixed_idx)

    U = np.zeros(ndof, dtype=float)
    if free_idx.size:
        rhs = F[free_idx]
        if fixed_idx.size:
            rhs -= K[np.ix_(free_idx, fixed_idx)] @ fixed_vals[fixed_idx]
        U[free_idx] = np.linalg.solve(K[np.ix_(free_idx, free_idx)], rhs)
    if fixed_idx.size:
        U[fixed_idx] = fixed_vals[fixed_idx]

    reactions = K @ U - F

    for ee, (n1, n2, n3) in enumerate(elements):
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]
        x3, y3 = nodes[n3]
        B, _ = _triangle_B_matrix(x1, y1, x2, y2, x3, y3)
        dofs = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n3, 2 * n3 + 1])
        ue = U[dofs]
        elem_sigma[ee, :] = D @ (B @ ue)

    return {
        "u": U,
        "disp": U.reshape(nn, 2),
        "reactions": reactions,
        "element_stress": elem_sigma,
        "nodes": nodes,
        "elements": elements,
    }
