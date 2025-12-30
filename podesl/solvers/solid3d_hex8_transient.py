from __future__ import annotations

from typing import Callable, Dict, Optional

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from .elastic3d import _D_iso_3d, _Ke_hex8


def _hex8_mass(xe: np.ndarray, rho: float) -> np.ndarray:
    gp = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
    Me = np.zeros((24, 24), dtype=float)
    for xi in gp:
        for eta in gp:
            for zeta in gp:
                N, dNdx, detJ = _hex8_shape_functions(xe, xi, eta, zeta)
                Nmat = np.zeros((3, 24))
                for a in range(8):
                    Nmat[0, 3 * a + 0] = N[a]
                    Nmat[1, 3 * a + 1] = N[a]
                    Nmat[2, 3 * a + 2] = N[a]
                Me += rho * (Nmat.T @ Nmat) * detJ
    return Me


def _hex8_shape_functions(xe: np.ndarray, xi: float, eta: float, zeta: float):
    # Reuse the same shape functions used in elastic3d._Ke_hex8.
    # The helper there is private; duplicating the minimal part here.
    N = np.array(
        [
            0.125 * (1 - xi) * (1 - eta) * (1 - zeta),
            0.125 * (1 + xi) * (1 - eta) * (1 - zeta),
            0.125 * (1 + xi) * (1 + eta) * (1 - zeta),
            0.125 * (1 - xi) * (1 + eta) * (1 - zeta),
            0.125 * (1 - xi) * (1 - eta) * (1 + zeta),
            0.125 * (1 + xi) * (1 - eta) * (1 + zeta),
            0.125 * (1 + xi) * (1 + eta) * (1 + zeta),
            0.125 * (1 - xi) * (1 + eta) * (1 + zeta),
        ],
        dtype=float,
    )
    dN = np.zeros((8, 3), dtype=float)
    dN[:, 0] = 0.125 * np.array(
        [
            -(1 - eta) * (1 - zeta),
            +(1 - eta) * (1 - zeta),
            +(1 + eta) * (1 - zeta),
            -(1 + eta) * (1 - zeta),
            -(1 - eta) * (1 + zeta),
            +(1 - eta) * (1 + zeta),
            +(1 + eta) * (1 + zeta),
            -(1 + eta) * (1 + zeta),
        ]
    )
    dN[:, 1] = 0.125 * np.array(
        [
            -(1 - xi) * (1 - zeta),
            -(1 + xi) * (1 - zeta),
            +(1 + xi) * (1 - zeta),
            +(1 - xi) * (1 - zeta),
            -(1 - xi) * (1 + zeta),
            -(1 + xi) * (1 + zeta),
            +(1 + xi) * (1 + zeta),
            +(1 - xi) * (1 + zeta),
        ]
    )
    dN[:, 2] = 0.125 * np.array(
        [
            -(1 - xi) * (1 - eta),
            -(1 + xi) * (1 - eta),
            -(1 + xi) * (1 + eta),
            -(1 - xi) * (1 + eta),
            +(1 - xi) * (1 - eta),
            +(1 + xi) * (1 - eta),
            +(1 + xi) * (1 + eta),
            +(1 - xi) * (1 + eta),
        ]
    )
    J = np.zeros((3, 3))
    for a in range(8):
        J[:, 0] += dN[a, 0] * xe[a, :]
        J[:, 1] += dN[a, 1] * xe[a, :]
        J[:, 2] += dN[a, 2] * xe[a, :]
    detJ = np.linalg.det(J)
    if detJ <= 0:
        raise ValueError("Invalid HEX8 element mapping.")
    invJ = np.linalg.inv(J)
    dNdx = dN @ invJ.T
    return N, dNdx, detJ


def _rayleigh_damping(K: csc_matrix, M: csc_matrix, ratio: float) -> csc_matrix:
    if ratio <= 0:
        return csc_matrix(K.shape)
    return ratio * (K + M)


def _assemble_hex8_KM(nodes: np.ndarray, elements: np.ndarray, E: float, nu: float, rho: float):
    nn = nodes.shape[0]
    ndof = 3 * nn
    K = lil_matrix((ndof, ndof), dtype=float)
    M = lil_matrix((ndof, ndof), dtype=float)
    D = _D_iso_3d(E, nu)
    for conn in elements:
        xe = nodes[conn, :]
        Ke, _ = _Ke_hex8(xe, D)
        Me = _hex8_mass(xe, rho)
        dofs = []
        for a in conn:
            dofs.extend([3 * a, 3 * a + 1, 3 * a + 2])
        for i, I in enumerate(dofs):
            for j, J in enumerate(dofs):
                K[I, J] += Ke[i, j]
                M[I, J] += Me[i, j]
    return csc_matrix(K), csc_matrix(M)


def solve_solid3d_hex8_transient(model: Dict[str, object]) -> Dict[str, np.ndarray]:
    nodes = np.asarray(model["nodes"], float)
    elements = np.asarray(model["elements"], int)
    E = float(model["E"])
    nu = float(model["nu"])
    rho = float(model["rho"])
    damping_ratio = float(model.get("damping_ratio", 0.0))
    fixed_dofs = np.asarray(model.get("fixed_dofs", []), dtype=int)
    load_func: Callable[[float], np.ndarray] = model.get("load_func", lambda t: np.zeros(3 * nodes.shape[0]))
    dt = float(model["dt"])
    t_end = float(model["t_end"])
    newmark = model.get("newmark", {}) or {}
    beta = float(newmark.get("beta", 0.25))
    gamma = float(newmark.get("gamma", 0.5))

    nn = nodes.shape[0]
    ndof = 3 * nn

    K, M = _assemble_hex8_KM(nodes, elements, E, nu, rho)
    C = _rayleigh_damping(K, M, damping_ratio)

    steps = int(round(t_end / dt))
    time = np.linspace(0.0, steps * dt, steps + 1)

    free = np.setdiff1d(np.arange(ndof), fixed_dofs)

    Kff = K[free][:, free]
    Mff = M[free][:, free]
    Cff = C[free][:, free]

    a0 = 1.0 / (beta * dt * dt)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2 * beta) - 1
    a4 = gamma / beta - 1
    a5 = dt * (gamma / (2 * beta) - 1)

    K_eff = (a0 * Mff + a1 * Cff + Kff).tocsc()
    solver = splu(K_eff)

    u = np.zeros((steps + 1, ndof))
    v = np.zeros((steps + 1, ndof))
    a = np.zeros((steps + 1, ndof))

    for i in range(steps):
        f_full = np.asarray(load_func(time[i + 1]), dtype=float).reshape(-1)
        f_free = f_full[free]
        rhs = (
            f_free
            + Mff @ (a0 * u[i, free] + a2 * v[i, free] + a3 * a[i, free])
            + Cff @ (a1 * u[i, free] + a4 * v[i, free] + a5 * a[i, free])
        )
        u_next_free = solver.solve(rhs)
        u[i + 1, free] = u_next_free
        v[i + 1, free] = a1 * (u[i + 1, free] - u[i, free]) - a4 * v[i, free] - a5 * a[i, free]
        a[i + 1, free] = a0 * (u[i + 1, free] - u[i, free]) - a2 * v[i, free] - a3 * a[i, free]

    return {"u": u, "v": v, "a": a, "time": time}
