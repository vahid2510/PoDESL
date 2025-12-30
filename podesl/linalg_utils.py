from __future__ import annotations

import numpy as np


def safe_solve(A: np.ndarray, b: np.ndarray, *, rcond: float = 1e-12) -> np.ndarray:
    """Solve Ax = b with a least-squares fallback for singular matrices."""

    A = np.array(A, dtype=float, copy=False)
    b = np.array(b, dtype=float, copy=False)
    try:
        return np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        sol, *_ = np.linalg.lstsq(A, b, rcond=rcond)
        return sol


def generalized_eig(
    A: np.ndarray,
    B: np.ndarray,
    *,
    sort: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the symmetric generalized eigenproblem A x = λ B x.

    Falls back to solving B^{-1} A when a Cholesky factorization is not available.
    """

    A = np.array(A, dtype=float, copy=True)
    B = np.array(B, dtype=float, copy=True)

    def _symmetrize(mat: np.ndarray) -> np.ndarray:
        return 0.5 * (mat + mat.T)

    flipped = False
    try:
        L = np.linalg.cholesky(B)
    except np.linalg.LinAlgError:
        try:
            L = np.linalg.cholesky(-B)
            A = -A
            B = -B
            flipped = True
        except np.linalg.LinAlgError:
            vals, vecs = np.linalg.eig(safe_solve(B, A))
            vals = np.real(vals)
            vecs = np.real(vecs)
            if sort:
                order = np.argsort(vals)
                vals = vals[order]
                vecs = vecs[:, order]
            return vals, vecs

    LinvA = np.linalg.solve(L, A)
    C = np.linalg.solve(L, LinvA.T).T  # L^{-1} A L^{-T}
    C = _symmetrize(C)
    vals, vecs = np.linalg.eigh(C)
    vecs = np.linalg.solve(L.T, vecs)
    vecs = np.real(vecs)
    vals = np.real(vals)
    if flipped:
        vals = -vals
    if sort:
        order = np.argsort(vals)
        vals = vals[order]
        vecs = vecs[:, order]
    return vals, vecs
