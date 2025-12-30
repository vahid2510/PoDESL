from __future__ import annotations

import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve


def _uniform_mesh(L: float, nel: int) -> np.ndarray:
    if nel <= 0:
        raise ValueError("nel must be positive.")
    return np.linspace(0.0, L, nel + 1, dtype=float)


def _equivalent_body_force(body_force: float, h: float) -> np.ndarray:
    return (body_force * h / 2.0) * np.array([1.0, 1.0], dtype=float)


def solve_bar1d_static(env):
    L = float(env.get("L"))
    A = float(env.get("A"))
    E = float(env.get("E"))
    nel = int(env.get("nel", 10))

    # Support both legacy and new parameter names.
    body_force = env.get("body_force", env.get("p", 0.0))
    body_force = 0.0 if body_force is None else float(body_force)
    traction_right = env.get("traction_right", env.get("right_force", 0.0))
    traction_right = 0.0 if traction_right is None else float(traction_right)

    left_bc = str(env.get("left", "fixed")).lower()
    right_bc = str(env.get("right", "fixed")).lower()

    left_disp = float(env.get("U_left", 0.0))
    right_disp = float(env.get("U_right", 0.0))

    point_loads = env.get("point_loads", []) or []
    if isinstance(point_loads, dict):
        point_loads = [point_loads]

    x = _uniform_mesh(L, nel)
    h = L / nel
    nn = x.size

    K = lil_matrix((nn, nn), dtype=float)
    f = np.zeros(nn, dtype=float)

    ke = (E * A / h) * np.array([[1.0, -1.0], [-1.0, 1.0]], dtype=float)
    fe_body = _equivalent_body_force(body_force, h)

    for e in range(nel):
        i, j = e, e + 1
        dofs = (i, j)
        for a in range(2):
            f[dofs[a]] += fe_body[a]
            for b in range(2):
                K[dofs[a], dofs[b]] += ke[a, b]

    for item in point_loads:
        if isinstance(item, (list, tuple)):
            P, xpos = float(item[1]), float(item[0])
        else:
            P = float(item.get("value", item.get("P", 0.0)))
            xpos = float(item.get("x", item.get("node", 0.0)))
        idx = int(round(xpos / h))
        idx = max(0, min(nn - 1, idx))
        f[idx] += P

    f[-1] += traction_right

    def apply_dirichlet(node: int, value: float) -> None:
        col = K[:, node].toarray().ravel()
        f[:] -= col * value
        K.rows[node] = [node]
        K.data[node] = [1.0]
        for r in range(nn):
            if r == node:
                continue
            row = K.rows[r]
            data = K.data[r]
            if node in row:
                idx = row.index(node)
                row.pop(idx)
                data.pop(idx)
        f[node] = value

    if left_bc in ("fixed", "dirichlet", "u"):
        apply_dirichlet(0, left_disp)
    if right_bc in ("fixed", "dirichlet", "u"):
        apply_dirichlet(nn - 1, right_disp)

    U = spsolve(csc_matrix(K), f)

    strain = np.zeros(nel, dtype=float)
    stress = np.zeros(nel, dtype=float)
    for e in range(nel):
        strain[e] = (U[e + 1] - U[e]) / h
        stress[e] = E * strain[e]

    return {"x": x, "u": U, "sigma": stress}
