from __future__ import annotations
import numpy as np

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

def solve_plane2d_q4_static(env):
    nodes = np.asarray(env.get('nodes'), float)
    elems = np.asarray(env.get('elems'), int)
    E = float(env.get('E'))
    nu = float(env.get('nu'))
    t = float(env.get('t', 1.0))
    plane = str(env.get('plane', 'stress')).lower()

    if plane.startswith('stress'):
        coef = E/(1.0-nu**2)
        D = coef*np.array([[1.0, nu, 0.0],
                           [nu, 1.0, 0.0],
                           [0.0, 0.0, (1.0-nu)/2.0]], float)
    else:
        coef = E/((1.0+nu)*(1.0-2.0*nu))
        D = coef*np.array([[1.0-nu, nu, 0.0],
                           [nu, 1.0-nu, 0.0],
                           [0.0, 0.0, (1.0-2.0*nu)/2.0]], float)

    nn = nodes.shape[0]
    ndof = 2*nn
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)

    gp = [-1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)]

    for conn in elems:
        xe = nodes[conn,:]
        Ke = np.zeros((8,8))
        for xi in gp:
            for eta in gp:
                N, dNdxi = _q4_shape(xi, eta)
                J = dNdxi.T @ xe
                detJ = np.linalg.det(J)
                if detJ <= 0:
                    raise ValueError("Invalid Q4 geometry")
                invJ = np.linalg.inv(J)
                dNdx = dNdxi @ invJ.T
                B = np.zeros((3,8))
                for i in range(4):
                    dNi_dx, dNi_dy = dNdx[i,0], dNdx[i,1]
                    B[0,2*i]   = dNi_dx
                    B[1,2*i+1] = dNi_dy
                    B[2,2*i]   = dNi_dy
                    B[2,2*i+1] = dNi_dx
                w = detJ * t
                Ke += B.T @ D @ B * w
        dofs = []
        for n in conn:
            dofs.extend([2*n, 2*n+1])
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J] += Ke[i,j]

    for ent in env.get('loads', []):
        n = int(ent[0]); Fx = float(ent[1]); Fy = float(ent[2])
        F[2*n]   += Fx
        F[2*n+1] += Fy

    fixed = {}
    for ent in env.get('fix', []):
        if len(ent) == 3:
            n = int(ent[0]); dof = str(ent[1]).lower(); val = float(ent[2])
            if dof in ('ux','both'):
                fixed[2*n] = val
            if dof in ('uy','both'):
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
        Ufree = np.linalg.solve(K[np.ix_(free, free)], rhs)
        U[free] = Ufree
    if fixed_idx.size > 0:
        U[fixed_idx] = np.array([fixed[i] for i in fixed_idx], float)

    R = K @ U - F
    return {'U': U, 'R': R, 'K': K}
