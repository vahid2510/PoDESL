from __future__ import annotations

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from ..linalg_utils import generalized_eig
from .bc_utils import FRAME3D_DOF_MAP, dirichlet_vector
from .input_utils import merged_dirichlet, merged_point_loads


def _local_axes(x1, y1, z1, x2, y2, z2):
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    L = float(np.sqrt(dx * dx + dy * dy + dz * dz))
    if L <= 0:
        raise ValueError("Zero-length frame element.")
    lx, ly, lz = dx / L, dy / L, dz / L
    vx = np.array([lx, ly, lz], dtype=float)
    ref = np.array([0.0, 0.0, 1.0], dtype=float)
    if np.isclose(abs(np.dot(vx, ref)), 1.0):
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    vy = ref - np.dot(ref, vx) * vx
    vy /= np.linalg.norm(vy)
    vz = np.cross(vx, vy)
    R = np.vstack([vx, vy, vz])
    return L, R


def _k_local(E, G, A, Iy, Iz, J, L):
    EA = E * A / L
    GJ = G * J / L
    EIy = E * Iy
    EIz = E * Iz

    k = np.zeros((12, 12), dtype=float)

    k[0, 0] = k[6, 6] = EA
    k[0, 6] = k[6, 0] = -EA

    c1 = 12 * EIz / (L ** 3)
    c2 = 6 * EIz / (L ** 2)
    c3 = 4 * EIz / L
    c4 = 2 * EIz / L

    k[1, 1] = k[7, 7] = c1
    k[1, 7] = k[7, 1] = -c1
    k[1, 5] = k[5, 1] = c2
    k[1, 11] = k[11, 1] = c2
    k[5, 5] = k[11, 11] = c3
    k[5, 7] = k[7, 5] = -c2
    k[5, 11] = k[11, 5] = c4
    k[7, 11] = k[11, 7] = -c2

    d1 = 12 * EIy / (L ** 3)
    d2 = 6 * EIy / (L ** 2)
    d3 = 4 * EIy / L
    d4 = 2 * EIy / L

    k[2, 2] = k[8, 8] = d1
    k[2, 8] = k[8, 2] = -d1
    k[2, 4] = k[4, 2] = -d2
    k[2, 10] = k[10, 2] = -d2
    k[4, 4] = k[10, 10] = d3
    k[4, 8] = k[8, 4] = d2
    k[4, 10] = k[10, 4] = d4
    k[8, 10] = k[10, 8] = d2

    k[3, 3] = k[9, 9] = GJ
    k[3, 9] = k[9, 3] = -GJ

    return k


def _m_local(rho, A, Iy, Iz, J, L):
    m = rho * A * L
    Irot = rho * (Iy + Iz) * L
    M = np.zeros((12, 12), dtype=float)
    for node in (0, 6):
        for d in range(3):
            M[node + d, node + d] = m / 2.0
    for node in (3, 9):
        M[node, node] = Irot / 2.0
    return M


def _assemble_frame3d(model):
    nodes = np.asarray(model["nodes"], float)
    elements = np.asarray(model["elements"], int)
    nn = nodes.shape[0]
    ndof = 6 * nn

    E = float(model["E"])
    G = float(model["G"])
    A = float(model["A"])
    Iy = float(model["Iy"])
    Iz = float(model["Iz"])
    J = float(model["J"])
    rho = float(model.get("rho", 0.0))

    K = lil_matrix((ndof, ndof), dtype=float)
    M = lil_matrix((ndof, ndof), dtype=float)

    for conn in elements:
        n1, n2 = conn
        x1, y1, z1 = nodes[n1]
        x2, y2, z2 = nodes[n2]
        L, R = _local_axes(x1, y1, z1, x2, y2, z2)
        k_loc = _k_local(E, G, A, Iy, Iz, J, L)
        m_loc = _m_local(rho, A, Iy, Iz, J, L)
        T = np.zeros((12, 12), dtype=float)
        for i in range(4):
            T[3 * i : 3 * i + 3, 3 * i : 3 * i + 3] = R
        k_glob = T.T @ k_loc @ T
        m_glob = T.T @ m_loc @ T
        dofs = np.array(
            [6 * n1, 6 * n1 + 1, 6 * n1 + 2, 6 * n1 + 3, 6 * n1 + 4, 6 * n1 + 5,
             6 * n2, 6 * n2 + 1, 6 * n2 + 2, 6 * n2 + 3, 6 * n2 + 4, 6 * n2 + 5]
        )
        for i, I in enumerate(dofs):
            for j, Jg in enumerate(dofs):
                K[I, Jg] += k_glob[i, j]
                M[I, Jg] += m_glob[i, j]
    return K, M


def solve_frame3d_beamcolumn_static(model: dict) -> dict:
    nodes = np.asarray(model["nodes"], float)
    K, _ = _assemble_frame3d(model)
    ndof = K.shape[0]
    load_entries = merged_point_loads(
        model.get("loads"),
        model.get("point_loads"),
        dof_order=("ux", "uy", "uz", "rx", "ry", "rz"),
    )
    F = np.zeros(ndof, dtype=float)
    for node, comp in load_entries:
        base = 6 * node
        for local, value in enumerate(comp):
            F[base + local] += value

    bc_entries = merged_dirichlet(
        model.get("fix"),
        model.get("bcs"),
        dof_map=FRAME3D_DOF_MAP,
    )
    fixed, vec = dirichlet_vector(
        bc_entries,
        ndof=ndof,
        ndof_per_node=6,
        dof_map=FRAME3D_DOF_MAP,
    )
    free = np.setdiff1d(np.arange(ndof), fixed)
    vals = vec[fixed] if fixed.size else np.zeros(0, dtype=float)

    U = np.zeros(ndof, dtype=float)
    if fixed.size:
        U[fixed] = vals

    Kff = csc_matrix(K[free][:, free])
    rhs = F[free]
    if fixed.size:
        rhs -= K[free][:, fixed] @ vals
    U[free] = splu(Kff).solve(rhs)

    reactions = K @ U - F
    return {"u": U, "reactions": reactions}


def solve_frame3d_beamcolumn_modal(model: dict) -> dict:
    nodes = np.asarray(model["nodes"], float)
    K, M = _assemble_frame3d(model)
    ndof = K.shape[0]
    nmodes = int(model.get("nmodes", 6))

    bc_entries = merged_dirichlet(
        model.get("fix"),
        model.get("bcs"),
        dof_map=FRAME3D_DOF_MAP,
    )
    fixed, _vec = dirichlet_vector(
        bc_entries,
        ndof=ndof,
        ndof_per_node=6,
        dof_map=FRAME3D_DOF_MAP,
    )
    free = np.setdiff1d(np.arange(ndof), fixed)
    Kff = K[free][:, free].toarray()
    Mff = M[free][:, free].toarray()

    vals, vecs = generalized_eig(Kff, Mff)
    order = np.argsort(vals)
    vals = vals[order][:nmodes]
    vecs = vecs[:, order][:, :nmodes]
    omega = np.sqrt(np.clip(vals, 0, None))
    freq = omega / (2 * np.pi)

    modes = np.zeros((ndof, nmodes))
    for i in range(nmodes):
        mode = np.zeros(ndof)
        mode[free] = vecs[:, i]
        modes[:, i] = mode

    return {"omega": omega, "freq": freq, "modes": modes}
