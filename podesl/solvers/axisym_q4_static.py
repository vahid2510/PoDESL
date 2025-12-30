from __future__ import annotations
import numpy as np

from ..linalg_utils import safe_solve
def _q4_shape(xi, eta):
    N = 0.25*np.array([(1-xi)*(1-eta),
                       (1+xi)*(1-eta),
                       (1+xi)*(1+eta),
                       (1-xi)*(1+eta)], float)
    dN = 0.25*np.array([[-(1-eta), -(1-xi)],
                        [ +(1-eta), -(1+xi)],
                        [ +(1+eta), +(1+xi)],
                        [ -(1+eta), +(1-xi)]], float)
    return N, dN

def solve_axisym_q4_static(env):
    nodes = np.asarray(env.get('nodes'), float)
    elems = np.asarray(env.get('elems'), int)
    E = float(env.get('E'))
    nu = float(env.get('nu'))

    coef = E/((1+nu)*(1-2*nu))
    D = coef*np.array([[1-nu,    nu,      nu,     0],
                       [nu,      1-nu,   nu,     0],
                       [nu,      nu,     1-nu,  0],
                       [0,       0,      0,     (1-2*nu)/2]], float)

    nn = nodes.shape[0]
    ndof = 2*nn
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)

    gp = [-1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)]

    for conn in elems:
        re = nodes[conn,0]
        ze = nodes[conn,1]
        Ke = np.zeros((8,8))
        for xi in gp:
            for eta in gp:
                N, dNdxi = _q4_shape(xi, eta)
                J = dNdxi.T @ np.column_stack((re, ze))
                detJ = np.linalg.det(J)
                if detJ <= 0:
                    raise ValueError("Invalid axisymmetric Q4 geometry")
                invJ = np.linalg.inv(J)
                dNdx = dNdxi @ invJ.T
                r_gp = float(N @ re)
                if r_gp <= 0:
                    raise ValueError("Non-positive radius in axisymmetric element")
                B = np.zeros((4, 8))
                for i in range(4):
                    dri = dNdx[i,0]; dzi = dNdx[i,1]; Ni = N[i]
                    B[0,2*i]   = dri
                    B[1,2*i+1] = dzi
                    B[2,2*i]   = Ni / r_gp
                    B[3,2*i]   = dzi
                    B[3,2*i+1] = dri
                w = detJ * 2*np.pi * r_gp
                Ke += B.T @ D @ B * w
        dofs = []
        for n in conn:
            dofs.extend([2*n, 2*n+1])
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J] += Ke[i,j]

    for ent in env.get('loads', []):
        n = int(ent[0]); Fr = float(ent[1]); Fz = float(ent[2])
        F[2*n]   += Fr
        F[2*n+1] += Fz

    fixed = {}
    for ent in env.get('fix', []):
        if len(ent) == 3:
            n = int(ent[0]); dof = str(ent[1]).lower(); val = float(ent[2])
            if dof in ("ur","both"):
                fixed[2*n] = val
            if dof in ("uz","both"):
                fixed[2*n+1] = val
        elif len(ent) == 2:
            n = int(ent[0]); val = float(ent[1])
            fixed[2*n] = val; fixed[2*n+1] = val

    fixed_idx = np.array(sorted(fixed.keys()), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed_idx)
    U = np.zeros(ndof)
    if free.size > 0:
        rhs = F[free]
        if fixed_idx.size > 0:
            rhs -= K[np.ix_(free, fixed_idx)] @ np.array([fixed[i] for i in fixed_idx], float)
        Ufree = safe_solve(K[np.ix_(free, free)], rhs)
        U[free] = Ufree
    if fixed_idx.size > 0:
        U[fixed_idx] = np.array([fixed[i] for i in fixed_idx], float)

    R = K @ U - F
    return {"U": U, "R": R, "K": K}
