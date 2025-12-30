from __future__ import annotations
import numpy as np
from ..linalg_utils import generalized_eig

def _truss2d_geom(x1,y1,x2,y2):
    L = float(np.hypot(x2-x1, y2-y1))
    if L <= 0:
        raise ValueError("Zero-length truss element")
    c = (x2-x1)/L
    s = (y2-y1)/L
    return c,s,L

def solve_truss2d_modal(env):
    nodes = np.asarray(env.get('nodes'), float)
    elems = np.asarray(env.get('elems'), int)
    E = float(env.get('E'))
    A = float(env.get('A'))
    rho = float(env.get('rho'))
    nmodes = int(env.get('nmodes', 6))

    nn = nodes.shape[0]
    ndof = 2*nn
    K = np.zeros((ndof, ndof))
    M = np.zeros((ndof, ndof))

    for e,(n1,n2) in enumerate(elems):
        x1,y1 = nodes[n1]
        x2,y2 = nodes[n2]
        c,s,L = _truss2d_geom(x1,y1,x2,y2)
        k_loc = (E*A/L) * np.array([[ 1,-1],
                                    [-1, 1]], float)
        T = np.array([[ c, s, 0, 0],
                      [ 0, 0, c, s]], float)
        k_gl = T.T @ k_loc @ T
        m = rho*A*L
        m_loc = (m/6.0)*np.array([[2,1],
                                  [1,2]], float)
        m_gl = T.T @ m_loc @ T
        dofs = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J] += k_gl[i,j]
                M[I,J] += m_gl[i,j]

    fixed = {}
    for ent in env.get('fix', []):
        if len(ent) == 3:
            n = int(ent[0]); dof = str(ent[1]).lower(); val = float(ent[2])
            if dof in ("ux","both"):
                fixed[2*n] = val
            if dof in ("uy","both"):
                fixed[2*n+1] = val
        elif len(ent) == 2:
            n = int(ent[0]); val = float(ent[1])
            fixed[2*n] = val; fixed[2*n+1] = val

    fixed_idx = np.array(sorted(fixed.keys()), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed_idx)

    Kf = K[np.ix_(free, free)]
    Mf = M[np.ix_(free, free)]

    w, vec = generalized_eig(Kf, Mf)

    omegas = np.sqrt(np.clip(w, 0.0, None))
    freqs = omegas/(2.0*np.pi)

    nm = min(nmodes, vec.shape[1])
    modes = np.zeros((ndof, nm))
    for k in range(nm):
        full = np.zeros(ndof)
        full[free] = vec[:,k]
        mk = full @ (M @ full)
        if mk > 0:
            full /= np.sqrt(mk)
        modes[:,k] = full

    return {"freq": freqs[:nm], "omega": omegas[:nm], "modes": modes, "K": K, "M": M}
