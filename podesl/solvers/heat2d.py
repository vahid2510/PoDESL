
from __future__ import annotations
import numpy as np

def _q4_shape(xi, eta):
    N = 0.25*np.array([(1-xi)*(1-eta),
                       (1+xi)*(1-eta),
                       (1+xi)*(1+eta),
                       (1-xi)*(1+eta)], dtype=float)
    dN_dxi = 0.25*np.array([[-(1-eta), -(1-xi)],
                            [ +(1-eta), -(1+xi)],
                            [ +(1+eta), +(1+xi)],
                            [ -(1+eta), +(1-xi)]], dtype=float)
    return N, dN_dxi

def _t3_B_A(xe):
    x1,y1 = xe[0]; x2,y2 = xe[1]; x3,y3 = xe[2]
    A = 0.5*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    if A <= 0: raise ValueError("Invalid/zero-area T3")
    b = np.array([y2-y3, y3-y1, y1-y2], float)
    c = np.array([x3-x2, x1-x3, x2-x1], float)
    B = (1.0/(2*A))*np.array([[b[0], b[1], b[2]],
                               [c[0], c[1], c[2]]])
    return B, A

def _edge_pairs(etype):
    return {0:(0,1),1:(1,2),2:(2,3),3:(3,0)} if etype=="Q4" else {0:(0,1),1:(1,2),2:(2,0)}

def solve_heat2d(env):
    nodes = np.asarray(env.get("nodes"), float)
    elems = np.asarray(env.get("elems"), int)
    etype = env.get("etype","Q4").upper()
    k = float(env.get("k"))
    q = float(env.get("q", 0.0))

    nn = nodes.shape[0]; ndof = nn
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    if etype=="Q4":
        gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
        for conn in elems:
            xe = nodes[conn,:]
            Ke = np.zeros((4,4)); fe = np.zeros(4)
            for xi in gp:
                for eta in gp:
                    N, dN_dxi = _q4_shape(xi, eta)
                    J = np.zeros((2,2))
                    for a in range(4):
                        J[0,0]+=dN_dxi[a,0]*xe[a,0]; J[0,1]+=dN_dxi[a,0]*xe[a,1]
                        J[1,0]+=dN_dxi[a,1]*xe[a,0]; J[1,1]+=dN_dxi[a,1]*xe[a,1]
                    detJ = np.linalg.det(J)
                    if detJ <= 0: raise ValueError("Invalid Q4 mapping (detJ<=0)")
                    invJ = np.linalg.inv(J)
                    dN_dx = dN_dxi @ invJ.T
                    B = dN_dx
                    Ke += k * (B @ B.T) * detJ
                    if q != 0.0:
                        fe += N * q * detJ
            for i, I in enumerate(conn):
                F[I] += fe[i]
                for j, Jj in enumerate(conn):
                    K[I,Jj] += Ke[i,j]
    elif etype=="T3":
        for conn in elems:
            xe = nodes[conn,:]
            B, A = _t3_B_A(xe)
            Ke = k * (B.T @ B) * A
            fe = (q * A / 3.0) * np.ones(3)
            for i, I in enumerate(conn):
                F[I] += fe[i]
                for j, Jj in enumerate(conn):
                    K[I,Jj] += Ke[i,j]
    else:
        raise ValueError("etype must be Q4 or T3")

    for spec in env.get("qedge", []):
        if isinstance(spec, dict):
            ee=int(spec["elem"]); edge=int(spec["edge"]); qn=float(spec["qn"])
        else:
            ee, edge, qn = int(spec[0]), int(spec[1]), float(spec[2])
        conn = elems[ee]; a,b = _edge_pairs(etype)[edge]
        n1, n2 = conn[a], conn[b]
        x1,y1 = nodes[n1]; x2,y2 = nodes[n2]
        L = float(np.hypot(x2-x1,y2-y1))
        F[n1] += -qn * L * 0.5
        F[n2] += -qn * L * 0.5

    for spec in env.get("hedge", []):
        ee, edge, h, Tinf = int(spec[0]), int(spec[1]), float(spec[2]), float(spec[3])
        conn = elems[ee]; a,b = _edge_pairs(etype)[edge]
        n1, n2 = conn[a], conn[b]
        x1,y1 = nodes[n1]; x2,y2 = nodes[n2]
        L = float(np.hypot(x2-x1,y2-y1))
        Ke = (h*L/6.0)*np.array([[2,1],[1,2]])
        Fe = (h*Tinf*L/2.0)*np.array([1.0,1.0])
        idx = [n1, n2]
        for i, I in enumerate(idx):
            F[I] += Fe[i]
            for j, Jj in enumerate(idx):
                K[I,Jj] += Ke[i,j]

    fixed = {}
    for n, val in env.get("fix", []):
        fixed[int(n)] = float(val)
    fixed_idx = np.array(sorted(fixed.keys()), dtype=int)
    free_idx = np.setdiff1d(np.arange(ndof), fixed_idx)

    T = np.zeros(ndof)
    if free_idx.size>0:
        rhs = F[free_idx].copy()
        if fixed_idx.size>0:
            rhs -= K[np.ix_(free_idx, fixed_idx)] @ np.array([fixed[i] for i in fixed_idx])
        T[free_idx] = np.linalg.solve(K[np.ix_(free_idx, free_idx)], rhs)
    for i in fixed_idx: T[i] = fixed[i]
    return {"T": T, "nodes": nodes, "elems": elems, "etype": etype}
