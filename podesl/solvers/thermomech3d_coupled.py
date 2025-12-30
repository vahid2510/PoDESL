from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from .bc_utils import SOLID3D_DOF_MAP, dirichlet_vector
from .heat3d_transient import solve_heat3d_transient
from .elastic3d import _D_iso_3d, _Ke_hex8, _hex8_B
from .input_utils import merged_dirichlet, merged_point_loads


def _gauss_points():
    a = 1.0 / np.sqrt(3.0)
    pts = [(-a, -a, -a), (a, -a, -a), (a, a, -a), (-a, a, -a), (-a, -a, a), (a, -a, a), (a, a, a), (-a, a, a)]
    return pts


def _assemble_hex8_stiffness(nodes: np.ndarray, elements: np.ndarray, E: float, nu: float) -> csc_matrix:
    ndof = 3 * nodes.shape[0]
    K = lil_matrix((ndof, ndof), dtype=float)
    D = _D_iso_3d(E, nu)
    for conn in elements:
        xe = nodes[conn, :]
        Ke, _ = _Ke_hex8(xe, D)
        dofs = []
        for n in conn:
            dofs.extend([3 * n, 3 * n + 1, 3 * n + 2])
        for i, I in enumerate(dofs):
            for j, J in enumerate(dofs):
                K[I, J] += Ke[i, j]
    return csc_matrix(K)


def _build_mechanical_loads(loads, point_loads, ndof: int) -> np.ndarray:
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


def _split_dirichlet(nn: int, *sources) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    merged = merged_dirichlet(*sources, dof_map=SOLID3D_DOF_MAP)
    idx, vec = dirichlet_vector(
        merged,
        ndof=3 * nn,
        ndof_per_node=3,
        dof_map=SOLID3D_DOF_MAP,
    )
    vals = vec[idx] if idx.size else np.zeros(0, dtype=float)
    free = np.setdiff1d(np.arange(3 * nn), idx)
    return free, idx, vals


def _thermal_force(nodes: np.ndarray, elements: np.ndarray, D: np.ndarray, alpha: float, dT_elem: np.ndarray) -> np.ndarray:
    ndof = 3 * nodes.shape[0]
    F = np.zeros(ndof, dtype=float)
    gauss = _gauss_points()
    for ee, conn in enumerate(elements):
        xe = nodes[conn, :]
        dofs = []
        for n in conn:
            dofs.extend([3 * n, 3 * n + 1, 3 * n + 2])
        fe = np.zeros(24, dtype=float)
        delta = dT_elem[ee]
        eps_th = alpha * delta * np.array([1, 1, 1, 0, 0, 0], dtype=float)
        for (xi, eta, zeta) in gauss:
            B, detJ = _hex8_B(xe, xi, eta, zeta)
            fe += B.T @ (D @ eps_th) * detJ
        for i, I in enumerate(dofs):
            F[I] += fe[i]
    return F


def solve_thermomechanical_3d_coupled(model: Dict[str, object]) -> Dict[str, object]:
    thermal_model = dict(model["thermal_model"])
    mechanical_model = model["mechanical_model"]
    alpha = float(model["alpha"])
    T_ref = float(model.get("T_ref", 0.0))
    dt = float(model["dt"])
    t_end = float(model["t_end"])

    # Ensure thermal solver stores the full history.
    thermal_model.setdefault("store_every", 1)
    thermal_res = solve_heat3d_transient(thermal_model)
    Thist = thermal_res.get("Thist")
    if Thist is None:
        Thist = np.tile(thermal_res["T"], (int(round(t_end / dt)) + 1, 1))

    mech_nodes = np.asarray(mechanical_model["nodes"], dtype=float)
    mech_elems = np.asarray(mechanical_model["elements"], dtype=int)
    E = float(mechanical_model["E"])
    nu = float(mechanical_model["nu"])

    K = _assemble_hex8_stiffness(mech_nodes, mech_elems, E, nu)
    F_ext = _build_mechanical_loads(
        mechanical_model.get("loads"),
        mechanical_model.get("point_loads"),
        K.shape[0],
    )
    free, fixed, vals = _split_dirichlet(
        mech_nodes.shape[0],
        mechanical_model.get("fix"),
        mechanical_model.get("bcs"),
    )
    Kff = K[free][:, free]
    solver = splu(Kff)
    D = _D_iso_3d(E, nu)

    timeline = thermal_res["times"] if "times" in thermal_res else np.linspace(0.0, t_end, Thist.shape[0])
    disp_hist = np.zeros((timeline.size, K.shape[0]), dtype=float)

    for step in range(timeline.size):
        nodal_T = Thist[step]
        elem_T = np.array([nodal_T[elem].mean() - T_ref for elem in mech_elems], dtype=float)
        F_th = _thermal_force(mech_nodes, mech_elems, D, alpha, elem_T)
        rhs = F_ext - F_th
        if fixed.size:
            rhs_free = rhs[free] - K[free][:, fixed] @ vals
        else:
            rhs_free = rhs[free]
        U = np.zeros(K.shape[0], dtype=float)
        if fixed.size:
            U[fixed] = vals
        U[free] = solver.solve(rhs_free)
        disp_hist[step, :] = U

    return {
        "temperature": Thist,
        "displacement": disp_hist,
        "time": timeline,
    }
