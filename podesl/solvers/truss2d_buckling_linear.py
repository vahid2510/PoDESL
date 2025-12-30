from __future__ import annotations
import numpy as np
from ..linalg_utils import safe_solve

def _truss2d_geom(x1,y1,x2,y2):
    L = float(np.hypot(x2-x1, y2-y1))
    if L <= 0:
        raise ValueError("Zero-length truss element")
    c = (x2-x1)/L
    s = (y2-y1)/L
    return c,s,L

def solve_truss2d_buckling_linear(env):
    nodes = np.asarray(env.get('nodes'), float)
    elems = np.asarray(env.get('elems'), int)
    E = float(env.get('E'))
    A = float(env.get('A'))
    Nref_in = env.get('Nref', 0.0)

    nn = nodes.shape[0]
    ndof = 2*nn
    K = np.zeros((ndof, ndof))
    Kg = np.zeros((ndof, ndof))

    if np.ndim(Nref_in) == 0:
        Nrefs = np.full(elems.shape[0], float(Nref_in))
    else:
        Nrefs = np.asarray(Nref_in, float).reshape(-1)
        if Nrefs.size != elems.shape[0]:
            raise ValueError("Nref length must equal number of elements")

    for e,(n1,n2) in enumerate(elems):
        x1,y1 = nodes[n1]
        x2,y2 = nodes[n2]
        c,s,L = _truss2d_geom(x1,y1,x2,y2)
        k_loc = (E*A/L) * np.array([[ 1,-1],
                                    [-1, 1]], float)
        T = np.array([[ c, s, 0, 0],
                      [ 0, 0, c, s]], float)
        k_gl = T.T @ k_loc @ T
        N = float(Nrefs[e])
        Kg_loc = (N/L) * np.array([[ 1,-1],
                                   [-1, 1]], float)
        Kg_gl = T.T @ Kg_loc @ T
        dofs = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J] += k_gl[i,j]
                Kg[I,J] += Kg_gl[i,j]

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
    Kg_f = Kg[np.ix_(free, free)]

    Af = safe_solve(Kg_f, Kf)
    lam, vec = np.linalg.eig(Af)
    lam = np.real(lam); vec = np.real(vec)
    order = np.argsort(lam)
    lam = lam[order]; vec = vec[:,order]

    nmodes = int(env.get("nmodes", 6))
    nm = min(nmodes, vec.shape[1])
    modes = np.zeros((ndof, nm))
    for k in range(nm):
        full = np.zeros(ndof)
        full[free] = vec[:,k]
        modes[:,k] = full

    return {"lambda": lam[:nm], "modes": modes, "K": K, "Kg": Kg}
