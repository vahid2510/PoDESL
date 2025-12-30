from __future__ import annotations

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from ..linalg_utils import generalized_eig
from .bc_utils import SOLID3D_DOF_MAP, dirichlet_vector
from .elastic3d import _D_iso_3d
from .input_utils import merged_dirichlet, merged_point_loads


def _hex20_shape(xi: float, eta: float, zeta: float):
    corner_signs = [
        (-1, -1, -1),
        (1, -1, -1),
        (1, 1, -1),
        (-1, 1, -1),
        (-1, -1, 1),
        (1, -1, 1),
        (1, 1, 1),
        (-1, 1, 1),
    ]

    edge_info = [
        ("xi", -1, -1),  # 8
        ("eta", 1, -1),
        ("xi", 1, -1),
        ("eta", -1, -1),
        ("zeta", -1, -1),
        ("zeta", 1, -1),
        ("zeta", 1, 1),
        ("zeta", -1, 1),
        ("xi", -1, 1),
        ("eta", 1, 1),
        ("xi", 1, 1),
        ("eta", -1, 1),
    ]

    N = np.zeros(20)
    dN = np.zeros((20, 3))

    for a, (sx, sy, sz) in enumerate(corner_signs):
        term = (1 + xi * sx) * (1 + eta * sy) * (1 + zeta * sz)
        phi = xi * sx + eta * sy + zeta * sz - 2
        N[a] = 0.125 * term * phi
        dphi = np.array([sx, sy, sz], dtype=float)
        dterm = np.array(
            [
                sx * (1 + eta * sy) * (1 + zeta * sz),
                sy * (1 + xi * sx) * (1 + zeta * sz),
                sz * (1 + xi * sx) * (1 + eta * sy),
            ],
            dtype=float,
        )
        dN[a, :] = 0.125 * (dterm * phi + term * dphi)

    for idx, info in enumerate(edge_info):
        a = 8 + idx
        if info[0] == "xi":
            sy, sz = info[1], info[2]
            N[a] = 0.25 * (1 - xi * xi) * (1 + eta * sy) * (1 + zeta * sz)
            dN[a, 0] = -0.5 * xi * (1 + eta * sy) * (1 + zeta * sz)
            dN[a, 1] = 0.25 * (1 - xi * xi) * sy * (1 + zeta * sz)
            dN[a, 2] = 0.25 * (1 - xi * xi) * (1 + eta * sy) * sz
        elif info[0] == "eta":
            sx, sz = info[1], info[2]
            N[a] = 0.25 * (1 + xi * sx) * (1 - eta * eta) * (1 + zeta * sz)
            dN[a, 0] = 0.25 * sx * (1 - eta * eta) * (1 + zeta * sz)
            dN[a, 1] = -0.5 * eta * (1 + xi * sx) * (1 + zeta * sz)
            dN[a, 2] = 0.25 * (1 + xi * sx) * (1 - eta * eta) * sz
        else:  # zeta edge
            sx, sy = info[1], info[2]
            N[a] = 0.25 * (1 + xi * sx) * (1 + eta * sy) * (1 - zeta * zeta)
            dN[a, 0] = 0.25 * sx * (1 + eta * sy) * (1 - zeta * zeta)
            dN[a, 1] = 0.25 * (1 + xi * sx) * sy * (1 - zeta * zeta)
            dN[a, 2] = -0.5 * zeta * (1 + xi * sx) * (1 + eta * sy)

    return N, dN


def _hex20_B(xe: np.ndarray, xi: float, eta: float, zeta: float):
    N, dN_dxi = _hex20_shape(xi, eta, zeta)
    J = np.zeros((3, 3))
    for a in range(20):
        J[:, 0] += dN_dxi[a, 0] * xe[a, :]
        J[:, 1] += dN_dxi[a, 1] * xe[a, :]
        J[:, 2] += dN_dxi[a, 2] * xe[a, :]
    detJ = np.linalg.det(J)
    if detJ <= 0:
        raise ValueError("Invalid HEX20 element geometry.")
    invJ = np.linalg.inv(J)
    dNdx = dN_dxi @ invJ.T
    B = np.zeros((6, 60), dtype=float)
    for a in range(20):
        dNx, dNy, dNz = dNdx[a, :]
        col = 3 * a
        B[0, col + 0] = dNx
        B[1, col + 1] = dNy
        B[2, col + 2] = dNz
        B[3, col + 0] = dNy
        B[3, col + 1] = dNx
        B[4, col + 1] = dNz
        B[4, col + 2] = dNy
        B[5, col + 0] = dNz
        B[5, col + 2] = dNx
    return N, B, detJ


def _gauss_3():
    a = np.sqrt(3.0 / 5.0)
    pts = [-a, 0.0, a]
    w = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
    return pts, w


def assemble_solid3d_hex20_KM(nodes: np.ndarray, elements: np.ndarray, E: float, nu: float, rho: float):
    nn = nodes.shape[0]
    ndof = 3 * nn
    K = lil_matrix((ndof, ndof), dtype=float)
    M = lil_matrix((ndof, ndof), dtype=float)
    D = _D_iso_3d(E, nu)
    gp, gw = _gauss_3()
    for conn in elements:
        xe = nodes[conn, :]
        dofs = []
        for n in conn:
            dofs.extend([3 * n, 3 * n + 1, 3 * n + 2])
        Ke = np.zeros((60, 60), dtype=float)
        Me = np.zeros((60, 60), dtype=float)
        for ix, wx in zip(gp, gw):
            for iy, wy in zip(gp, gw):
                for iz, wz in zip(gp, gw):
                    N, B, detJ = _hex20_B(xe, ix, iy, iz)
                    weight = wx * wy * wz * detJ
                    Ke += B.T @ D @ B * weight
                    Nmat = np.zeros((3, 60), dtype=float)
                    for a in range(20):
                        Nmat[0, 3 * a + 0] = N[a]
                        Nmat[1, 3 * a + 1] = N[a]
                        Nmat[2, 3 * a + 2] = N[a]
                    Me += rho * (Nmat.T @ Nmat) * weight
        for i, I in enumerate(dofs):
            for j, J in enumerate(dofs):
                K[I, J] += Ke[i, j]
                M[I, J] += Me[i, j]
    return csc_matrix(K), csc_matrix(M)


def _build_load_vector(loads, point_loads, ndof):
    entries = merged_point_loads(
        loads,
        point_loads,
        dof_order=("ux", "uy", "uz"),
    )
    F = np.zeros(ndof, dtype=float)
    for node, comp in entries:
        base = 3 * node
        for local, value in enumerate(comp[:3]):
            F[base + local] += value
    return F


def solve_solid3d_hex20_static(model: dict) -> dict:
    nodes = np.asarray(model["nodes"], float)
    elements = np.asarray(model["elements"], int)
    E = float(model["E"])
    nu = float(model["nu"])
    rho = float(model.get("rho", 0.0))
    ndof = nodes.shape[0] * 3

    K, _ = assemble_solid3d_hex20_KM(nodes, elements, E, nu, rho)
    F = _build_load_vector(model.get("loads"), model.get("point_loads"), ndof)
    bc_entries = merged_dirichlet(
        model.get("fix"),
        model.get("bcs"),
        dof_map=SOLID3D_DOF_MAP,
    )
    fixed, vec = dirichlet_vector(
        bc_entries,
        ndof=ndof,
        ndof_per_node=3,
        dof_map=SOLID3D_DOF_MAP,
    )
    free = np.setdiff1d(np.arange(ndof), fixed)
    vals = vec[fixed] if fixed.size else np.zeros(0, dtype=float)

    U = np.zeros(ndof, dtype=float)
    if fixed.size:
        U[fixed] = vals

    K_ff = K[free][:, free]
    rhs = F[free]
    if fixed.size:
        rhs -= K[free][:, fixed] @ vals
    U[free] = splu(K_ff).solve(rhs)
    reactions = K @ U - F
    return {"u": U, "reactions": reactions}


def solve_solid3d_hex20_modal(model: dict) -> dict:
    nodes = np.asarray(model["nodes"], float)
    elements = np.asarray(model["elements"], int)
    E = float(model["E"])
    nu = float(model["nu"])
    rho = float(model["rho"])
    ndof = nodes.shape[0] * 3
    nmodes = int(model.get("nmodes", 6))

    K, M = assemble_solid3d_hex20_KM(nodes, elements, E, nu, rho)
    bc_entries = merged_dirichlet(
        model.get("fix"),
        model.get("bcs"),
        dof_map=SOLID3D_DOF_MAP,
    )
    fixed, _vec = dirichlet_vector(
        bc_entries,
        ndof=ndof,
        ndof_per_node=3,
        dof_map=SOLID3D_DOF_MAP,
    )
    free = np.setdiff1d(np.arange(ndof), fixed)
    K_ff = K[free][:, free]
    M_ff = M[free][:, free]

    vals, vecs = generalized_eig(K_ff.toarray(), M_ff.toarray())
    idx = np.argsort(vals)[:nmodes]
    vals = vals[idx]
    vecs = vecs[:, idx]
    omega = np.sqrt(np.clip(vals, 0, None))
    freq = omega / (2 * np.pi)

    modes = np.zeros((ndof, len(idx)))
    for i in range(len(idx)):
        modes[free, i] = vecs[:, i]

    return {"omega": omega, "freq": freq, "modes": modes}
