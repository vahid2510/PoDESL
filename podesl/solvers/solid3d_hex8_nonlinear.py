from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu

from .elastic3d import _D_iso_3d, _hex8_B
from .bc_utils import SOLID3D_DOF_MAP, dirichlet_vector
from .input_utils import merged_dirichlet, merged_point_loads


def _gauss_points():
    a = 1.0 / np.sqrt(3.0)
    pts = [(-a, -a, -a), (a, -a, -a), (a, a, -a), (-a, a, -a), (-a, -a, a), (a, -a, a), (a, a, a), (-a, a, a)]
    return pts


@dataclass
class NonlinearMaterial:
    E: float
    nu: float
    k: float

    @property
    def D_linear(self) -> np.ndarray:
        return _D_iso_3d(self.E, self.nu)

    def stress_tangent(self, eps: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        D = self.D_linear
        if self.k == 0.0:
            return D @ eps, D
        nonlinear = self.k * np.power(eps, 3)
        tangent_extra = np.diag(3.0 * self.k * np.power(eps, 2))
        return D @ eps + nonlinear, D + tangent_extra


def _build_load_vector(loads, point_loads, ndof: int) -> np.ndarray:
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


def _assemble_tangent_residual(nodes: np.ndarray, elements: np.ndarray, U: np.ndarray, material: NonlinearMaterial):
    nn = nodes.shape[0]
    ndof = 3 * nn
    K = lil_matrix((ndof, ndof), dtype=float)
    Fint = np.zeros(ndof, dtype=float)

    gauss = _gauss_points()
    for conn in elements:
        xe = nodes[conn, :]
        dofs = []
        for n in conn:
            dofs.extend([3 * n, 3 * n + 1, 3 * n + 2])
        ue = U[dofs]
        Ke = np.zeros((24, 24), dtype=float)
        fint = np.zeros(24, dtype=float)

        for (xi, eta, zeta) in gauss:
            B, detJ = _hex8_B(xe, xi, eta, zeta)
            eps = B @ ue
            sigma, D_tan = material.stress_tangent(eps)
            Ke += B.T @ D_tan @ B * detJ
            fint += B.T @ sigma * detJ

        for i, I in enumerate(dofs):
            Fint[I] += fint[i]
            for j, J in enumerate(dofs):
                K[I, J] += Ke[i, j]
    return K, Fint


def solve_solid3d_hex8_nonlinear(model: Dict[str, object]) -> Dict[str, object]:
    nodes = np.asarray(model["nodes"], float)
    elements = np.asarray(model["elements"], int)
    material_data = model.get("material", {}) or {}
    material = NonlinearMaterial(
        E=float(material_data.get("E", model.get("E", 210e9))),
        nu=float(material_data.get("nu", model.get("nu", 0.3))),
        k=float(material_data.get("nonlinear_coef", 0.0)),
    )

    loads = model.get("loads")
    point_loads = model.get("point_loads")
    bcs = model.get("bcs")
    max_iter = int(model.get("max_iter", 25))
    tol = float(model.get("tol", 1e-6))
    load_steps = int(model.get("load_steps", 10))

    nn = nodes.shape[0]
    ndof = 3 * nn

    F_total = _build_load_vector(loads, point_loads, ndof)

    bc_entries = merged_dirichlet(
        bcs,
        model.get("fix"),
        dof_map=SOLID3D_DOF_MAP,
    )
    bc_dofs, bc_vec = dirichlet_vector(
        bc_entries,
        ndof=ndof,
        ndof_per_node=3,
        dof_map=SOLID3D_DOF_MAP,
    )
    bc_vals = bc_vec[bc_dofs] if bc_dofs.size else np.zeros(0, dtype=float)
    free = np.setdiff1d(np.arange(ndof), bc_dofs)

    U = np.zeros(ndof, dtype=float)
    if bc_dofs.size:
        U[bc_dofs] = bc_vals
    u_steps: List[np.ndarray] = []
    iterations: List[int] = []
    converged = True

    res_history = []
    line_search = bool(model.get("line_search", False))
    for step in range(1, load_steps + 1):
        lam = step / load_steps
        F_target = lam * F_total
        for it in range(1, max_iter + 1):
            K, Fint = _assemble_tangent_residual(nodes, elements, U, material)
            R = F_target - Fint
            R_free = R[free]
            res_norm = np.linalg.norm(R_free, ord=np.inf)
            res_history.append(res_norm)
            if res_norm < tol:
                iterations.append(it)
                break
            K_ff = csc_matrix(K[free][:, free])
            delta = splu(K_ff).solve(R_free)
            alpha = 1.0
            if line_search:
                for _ in range(5):
                    trial = U.copy()
                    trial[free] += alpha * delta
                    _, Fint_trial = _assemble_tangent_residual(nodes, elements, trial, material)
                    R_trial = F_target - Fint_trial
                    if np.linalg.norm(R_trial[free], ord=np.inf) < res_norm:
                        break
                    alpha *= 0.5
            U[free] += alpha * delta
            if bc_dofs.size:
                U[bc_dofs] = bc_vals
            if it == max_iter:
                iterations.append(it)
                converged = False
        u_steps.append(U.copy())
        if not converged:
            break

    return {
        "u": U,
        "u_steps": u_steps,
        "iterations": iterations,
        "residual_history": res_history,
        "converged": converged,
    }
