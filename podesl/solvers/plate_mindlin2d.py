
from __future__ import annotations
import numpy as np

def _q4_shape(xi, eta):
    N = 0.25*np.array([
        (1-xi)*(1-eta),
        (1+xi)*(1-eta),
        (1+xi)*(1+eta),
        (1-xi)*(1+eta),
    ])
    dN_dxi = 0.25*np.array([
        [-(1-eta), -(1-xi)],
        [ +(1-eta), -(1+xi)],
        [ +(1+eta), +(1+xi)],
        [ -(1+eta), +(1-xi)],
    ])
    return N, dN_dxi

def _jacobian_Q4(xe, dN_dxi):
    J = np.zeros((2,2))
    for a in range(4):
        J[0,0] += dN_dxi[a,0]*xe[a,0]; J[0,1] += dN_dxi[a,0]*xe[a,1]
        J[1,0] += dN_dxi[a,1]*xe[a,0]; J[1,1] += dN_dxi[a,1]*xe[a,1]
    detJ = np.linalg.det(J)
    if detJ <= 0:
        raise ValueError("Invalid Q4 mapping (detJ<=0). Check node order.")
    invJ = np.linalg.inv(J)
    dN_dx = np.zeros((4,2))
    for a in range(4):
        dN_dx[a,:] = invJ @ dN_dxi[a,:]
    return dN_dx, detJ

def solve_plate_mindlin_q4(env):
    nodes = np.asarray(env.get("nodes"), dtype=float)
    elems = np.asarray(env.get("elems"), dtype=int)
    E = float(env.get("E", 210e9))
    nu = float(env.get("nu", 0.3))
    t  = float(env.get("t", 0.01))
    kappa = float(env.get("kappa", 5.0/6.0))
    q = float(env.get("q", 0.0))

    nn = nodes.shape[0]
    ndof = 3*nn

    D = (E * t**3) / (12.0 * (1.0 - nu**2))
    Db = D * np.array([
        [1.0,    nu,        0.0],
        [nu,     1.0,       0.0],
        [0.0,    0.0, (1.0-nu)/2.0]
    ])
    G = E/(2.0*(1.0+nu))
    Ds = kappa * G * t * np.eye(2)

    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)

    gp2 = [-1/np.sqrt(3), 1/np.sqrt(3)]
    gp1 = [0.0]

    for e, conn in enumerate(elems):
        xe = nodes[conn,:]
        Ke = np.zeros((12,12))
        Fe = np.zeros(12)

        for xi in gp2:
            for eta in gp2:
                N, dN_dxi = _q4_shape(xi, eta)
                dN_dx, detJ = _jacobian_Q4(xe, dN_dxi)
                Bb = np.zeros((3,12))
                for a in range(4):
                    Bb[0, 3*a+1] = dN_dx[a,0]
                    Bb[1, 3*a+2] = dN_dx[a,1]
                    Bb[2, 3*a+1] = dN_dx[a,1]
                    Bb[2, 3*a+2] = dN_dx[a,0]
                Ke += (Bb.T @ Db @ Bb) * detJ

        for xi in gp1:
            for eta in gp1:
                N, dN_dxi = _q4_shape(xi, eta)
                dN_dx, detJ = _jacobian_Q4(xe, dN_dxi)
                Bs = np.zeros((2,12))
                for a in range(4):
                    Bs[0, 3*a + 0] = dN_dx[a,0]
                    Bs[1, 3*a + 0] = dN_dx[a,1]
                    Bs[0, 3*a + 1] = N[a]
                    Bs[1, 3*a + 2] = N[a]
                Ke += (Bs.T @ Ds @ Bs) * detJ
                if q != 0.0:
                    Fe_w = (q * detJ) * np.ones(4) / 4.0
                    for a in range(4):
                        Fe[3*a + 0] += Fe_w[a]

        dofs = []
        for a in conn:
            dofs += [3*a+0, 3*a+1, 3*a+2]
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J] += Ke[i,j]
        for i,I in enumerate(dofs):
            F[I] += Fe[i]

    for (node, pz) in env.get("point_loads", []):
        F[3*node + 0] += float(pz)

    fixed = {}
    for (node, val) in env.get("fixW", []):
        fixed[3*int(node)+0] = float(val)
    for (node, comp, val) in env.get("fixRot", []):
        node = int(node); val = float(val)
        if comp == 'rx': fixed[3*node+1] = val
        elif comp == 'ry': fixed[3*node+2] = val
        else: raise ValueError("fixRot comp must be 'rx' or 'ry'")

    all_idx = np.arange(ndof)
    fixed_idx = np.array(sorted(fixed.keys()), dtype=int)
    free_idx = np.setdiff1d(all_idx, fixed_idx)

    U = np.zeros(ndof)
    if fixed_idx.size > 0:
        U[fixed_idx] = np.array([fixed[i] for i in fixed_idx])
        F[free_idx] -= K[np.ix_(free_idx, fixed_idx)] @ U[fixed_idx]
        Kff = K[np.ix_(free_idx, free_idx)]
        Ff = F[free_idx]
        if Kff.size > 0:
            U[free_idx] = np.linalg.solve(Kff, Ff)
    else:
        U[:] = np.linalg.solve(K, F)

    wnod = U[0::3]
    rx = U[1::3]
    ry = U[2::3]

    M_elem = []
    for conn in elems:
        xe = nodes[conn,:]
        N, dN_dxi = _q4_shape(0.0, 0.0)
        dN_dx, detJ = _jacobian_Q4(xe, dN_dxi)
        kx = 0.0; ky = 0.0; kxy = 0.0
        for i,a in enumerate(conn):
            kx  += dN_dx[a,0] * rx[a]
            ky  += dN_dx[a,1] * ry[a]
            kxy += dN_dx[a,1] * rx[a] + dN_dx[a,0] * ry[a]
        D = (E * t**3) / (12.0 * (1.0 - nu**2))
        Db = D * np.array([[1.0, nu, 0.0],[nu, 1.0, 0.0],[0.0, 0.0, (1.0-nu)/2.0]])
        M = Db @ np.array([kx,ky,kxy])
        M_elem.append(M)
    M_elem = np.asarray(M_elem)

    return {"w": wnod, "theta": np.vstack([rx,ry]).T, "nodes": nodes, "elems": elems, "M_elem": M_elem}
