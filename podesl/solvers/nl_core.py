from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple

try:
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    _HAS_SPARSE = True
except Exception:
    sp = None; spla = None; _HAS_SPARSE = False

from ..linalg_utils import safe_solve

@dataclass
class NLControls:
    method: str = "loadstep"  # "loadstep" | "arc"
    load_steps: int = 20
    max_iter: int = 40
    tol: float = 1e-9
    line_search: bool = True
    arc_len: float = 1e-3
    arc_minmax: Tuple[float,float] = (1e-5, 1e-1)
    report_each: int = 1
    sparse: bool = False

def _solve_linear(A, b, free, controls: NLControls):
    if controls.sparse and _HAS_SPARSE and sp.issparse(A):
        return spla.spsolve(A[free][:,free], b[free])
    if controls.sparse and _HAS_SPARSE and not sp.issparse(A):
        A = sp.csr_matrix(A)
        return spla.spsolve(A[free][:,free], b[free])
    return safe_solve(A[free][:,free], b[free])

def _backtracking(obj_fn, alpha0=1.0, c=1e-4, rho=0.5, max_back=10):
    alpha = alpha0
    f0 = obj_fn(0.0)
    for _ in range(max_back):
        f1 = obj_fn(alpha)
        if f1 <= (1.0 - c*alpha)*f0:
            return alpha
        alpha *= rho
    return 0.0

def newton_loadstep(assemble_cb, U0, Fext, free, controls: NLControls):
    U = U0.copy()
    history = []
    for it in range(controls.max_iter):
        K, Fint = assemble_cb(U)
        R = Fext - Fint
        res = np.linalg.norm(R[free], ord=np.inf)
        history.append(float(res))
        if res < controls.tol: return U, {"iters": it+1, "res_hist": history, "status":"converged"}
        du_free = _solve_linear(K, R, free, controls)
        if controls.line_search:
            def phi(a):
                Ua = U.copy(); Ua[free] += a*du_free
                _, Fint_a = assemble_cb(Ua)
                r = Fext - Fint_a
                return float(0.5*np.dot(r[free], r[free]))
            a = _backtracking(phi, alpha0=1.0)
        else:
            a = 1.0
        U[free] += a*du_free
        if a == 0.0: break
    return U, {"iters": controls.max_iter, "res_hist": history, "status":"stalled"}

def solve_load_steps(assemble_cb, ndof, free, Fext_full, controls: NLControls):
    U = np.zeros(ndof)
    path = []
    lam = 0.0
    dlam = 1.0/controls.load_steps
    for step in range(1, controls.load_steps+1):
        lam_t = min(1.0, lam + dlam)
        Fext = lam_t * Fext_full
        U, info = newton_loadstep(assemble_cb, U, Fext, free, controls)
        lam = lam_t
        path.append([step, lam, info["res_hist"][-1] if info["res_hist"] else 0.0])
    return U, np.array(path), {"status":"ok"}

def solve_arc_length(assemble_cb, ndof, free, Fext_full, controls: NLControls):
    U = np.zeros(ndof); lam = 0.0
    s = controls.arc_len
    dU_prev = np.zeros(ndof)
    path = []
    for step in range(1, controls.load_steps+1):
        for it in range(controls.max_iter):
            K, Fint = assemble_cb(U)
            R = lam*Fext_full - Fint
            F = Fext_full
            A = np.zeros((K.shape[0]+1, K.shape[1]+1))
            A[:ndof,:ndof] = K
            A[:ndof,-1] = -F
            v = dU_prev
            A[-1,:ndof] = 2*v
            A[-1,-1] = 0.0
            rhs_aug = np.zeros(ndof+1)
            rhs_aug[:ndof] = -R
            rhs_aug[-1] = s*s - np.dot(v,v)
            sol = safe_solve(A, rhs_aug)
            dU = sol[:ndof]; dlam = sol[-1]
            def phi(a):
                Ua = U + a*dU; lama = lam + a*dlam
                _, Fint_a = assemble_cb(Ua)
                r = lama*Fext_full - Fint_a
                return float(0.5*np.dot(r[free], r[free]))
            a = _backtracking(phi, alpha0=1.0)
            U = U + a*dU; lam = lam + a*dlam
            _, Fint2 = assemble_cb(U)
            r2 = lam*Fext_full - Fint2
            res = np.linalg.norm(r2[free], ord=np.inf)
            if res < controls.tol:
                dU_prev = a*dU
                if it+1 <= 4: s = min(controls.arc_minmax[1], s*1.5)
                elif it+1 >= 12: s = max(controls.arc_minmax[0], s*0.5)
                break
        path.append([step, lam, res])
        if lam >= 1.0 - 1e-12: break
    return U, np.array(path), {"status":"ok"}
