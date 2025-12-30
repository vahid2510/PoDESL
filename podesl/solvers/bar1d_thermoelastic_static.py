from __future__ import annotations
import numpy as np

def solve_bar1d_thermoelastic_static(env):
    nodes = np.asarray(env.get('nodes'), float).reshape(-1)
    elems = np.asarray(env.get('elems'), int)
    E = float(env.get('E'))
    A = float(env.get('A'))
    alpha = float(env.get('alpha'))
    T = np.asarray(env.get('T'), float).reshape(-1)
    Tref = float(env.get('Tref', 0.0))

    nn = nodes.shape[0]
    ndof = nn
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)
    F_th = np.zeros(ndof)

    for e,(n1,n2) in enumerate(elems):
        x1 = nodes[n1]; x2 = nodes[n2]
        L = float(x2 - x1)
        if L <= 0:
            raise ValueError("Non-positive element length in bar1d_thermoelastic")
        ke = (E*A/L) * np.array([[1.0,-1.0],
                                 [-1.0,1.0]], float)
        dofs = [n1, n2]
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J] += ke[i,j]

        Tavg = 0.5*(T[n1] + T[n2])
        dT = Tavg - Tref
        eps_th = alpha * dT
        fe_th = E*A*eps_th * np.array([-1.0, 1.0], float)
        for i,I in enumerate(dofs):
            F_th[I] += fe_th[i]

    for ent in env.get('loads', []):
        n = int(ent[0]); Fx = float(ent[1])
        F[n] += Fx

    fixed = {}
    for ent in env.get('fix', []):
        if len(ent) == 2:
            n = int(ent[0]); val = float(ent[1])
            fixed[n] = val
        elif len(ent) == 3:
            n = int(ent[0]); val = float(ent[2])
            fixed[n] = val

    fixed_idx = np.array(sorted(fixed.keys()), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed_idx)

    U = np.zeros(ndof)
    Feff = F - F_th
    if free.size > 0:
        rhs = Feff[free]
        if fixed_idx.size > 0:
            rhs -= K[np.ix_(free, fixed_idx)] @ np.array([fixed[i] for i in fixed_idx], float)
        Ufree = np.linalg.solve(K[np.ix_(free, free)], rhs)
        U[free] = Ufree
    if fixed_idx.size > 0:
        U[fixed_idx] = np.array([fixed[i] for i in fixed_idx], float)

    R = K @ U - Feff
    return {"U": U, "R": R, "K": K}
