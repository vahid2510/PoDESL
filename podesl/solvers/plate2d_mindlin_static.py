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

def solve_plate2d_mindlin_static(env):
    nodes = np.asarray(env.get('nodes'), float)
    elems = np.asarray(env.get('elems'), int)
    E = float(env.get('E'))
    nu = float(env.get('nu'))
    t = float(env.get('t', 1.0))
    kappa = float(env.get('kappa', 5.0/6.0))

    nn = nodes.shape[0]
    ndof = 3*nn
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)

    D0 = E * t**3 / (12.0*(1.0-nu**2))
    Db = D0*np.array([[1.0, nu, 0.0],
                      [nu, 1.0, 0.0],
                      [0.0, 0.0, (1.0-nu)/2.0]], float)

    G = E / (2.0*(1.0+nu))
    Ds = kappa * G * t * np.eye(2)

    gp = [-1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)]

    for conn in elems:
        xe = nodes[conn,:]
        Ke = np.zeros((12,12))
        for xi in gp:
            for eta in gp:
                N, dNdxi = _q4_shape(xi, eta)
                J = dNdxi.T @ xe
                detJ = np.linalg.det(J)
                if detJ <= 0:
                    raise ValueError("Invalid Mindlin Q4 geometry")
                invJ = np.linalg.inv(J)
                dNdx = dNdxi @ invJ.T

                Bb = np.zeros((3,12))
                Bs = np.zeros((2,12))
                for i in range(4):
                    dNi_dx, dNi_dy = dNdx[i,0], dNdx[i,1]
                    Ni = N[i]
                    col_w  = 3*i
                    col_tx = 3*i+1
                    col_ty = 3*i+2

                    Bb[0, col_tx] = dNi_dx
                    Bb[1, col_ty] = dNi_dy
                    Bb[2, col_tx] = dNi_dy
                    Bb[2, col_ty] = dNi_dx

                    Bs[0, col_w]  = dNi_dx
                    Bs[0, col_tx] = Ni
                    Bs[1, col_w]  = dNi_dy
                    Bs[1, col_ty] = Ni
                wgt = detJ
                Ke += Bb.T @ Db @ Bb * wgt + Bs.T @ Ds @ Bs * wgt

        dofs = []
        for n in conn:
            dofs.extend([3*n, 3*n+1, 3*n+2])
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J] += Ke[i,j]

    for ent in env.get('loads', []):
        n = int(ent[0]); Fw = float(ent[1])
        F[3*n] += Fw

    fixed = {}
    for ent in env.get('fix', []):
        if len(ent) == 3:
            n = int(ent[0]); dof = str(ent[1]).lower(); val = float(ent[2])
            if dof in ('w','all','both'):
                fixed[3*n] = val
            if dof in ('tx','all'):
                fixed[3*n+1] = val
            if dof in ('ty','all'):
                fixed[3*n+2] = val
        elif len(ent) == 2:
            n = int(ent[0]); val = float(ent[1])
            fixed[3*n] = val; fixed[3*n+1] = val; fixed[3*n+2] = val

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
