from __future__ import annotations
import numpy as np

# 2D scalar heat conduction, Q4 element, steady state
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

def solve_heat2d_q4_static(env):
    nodes = np.asarray(env.get('nodes'), float)
    elems = np.asarray(env.get('elems'), int)
    k = float(env.get('k'))
    q = env.get('q', None)

    nn = nodes.shape[0]; ndof = nn
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    gp = [-1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)]
    for conn in elems:
        xe = nodes[conn,:]  # (4,2)
        Ke = np.zeros((4,4)); Fe = np.zeros(4)
        for xi in gp:
            for eta in gp:
                N, dNdxi = _q4_shape(xi, eta)
                J = dNdxi.T @ xe
                detJ = np.linalg.det(J)
                if detJ <= 0:
                    raise ValueError('Invalid Q4 geometry')
                invJ = np.linalg.inv(J)
                dNdx = dNdxi @ invJ.T   # (4,2)
                B = dNdx
                w = detJ
                Ke += k * (B @ B.T) * w
                if q is not None:
                    qv = float(q)
                    Fe += qv * N * w
        dofs = [int(n) for n in conn]
        for i,I in enumerate(dofs):
            F[I] += Fe[i]
            for j,J in enumerate(dofs):
                K[I,J] += Ke[i,j]

    # Neumann flux on edges: qedge = [[elem, edge, qn], ...]
    # edges: 0-1, 1-2, 2-3, 3-0
    for spec in env.get('qedge', []):
        ee = int(spec[0]); edge = int(spec[1]); qn = float(spec[2])
        conn = elems[ee]
        a,b = [(0,1),(1,2),(2,3),(3,0)][edge]
        n1, n2 = conn[a], conn[b]
        x1,y1 = nodes[n1]; x2,y2 = nodes[n2]
        L = float(np.hypot(x2-x1, y2-y1))
        Fe = qn * L/2.0 * np.array([1.0, 1.0], float)
        F[n1] += Fe[0]; F[n2] += Fe[1]

    fixed = {}
    for ent in env.get('fixT', []):
        if len(ent) == 2:
            n = int(ent[0]); val = float(ent[1]); fixed[n] = val
    fixed_idx = np.array(sorted(fixed.keys()), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed_idx)
    T = np.zeros(ndof)
    if free.size > 0:
        rhs = F[free]
        if fixed_idx.size > 0:
            rhs -= K[np.ix_(free, fixed_idx)] @ np.array([fixed[i] for i in fixed_idx], float)
        Tfree = np.linalg.solve(K[np.ix_(free, free)], rhs)
        T[free] = Tfree
    if fixed_idx.size > 0:
        T[fixed_idx] = np.array([fixed[i] for i in fixed_idx], float)

    R = K @ T - F
    return {'T': T, 'R': R, 'K': K}
