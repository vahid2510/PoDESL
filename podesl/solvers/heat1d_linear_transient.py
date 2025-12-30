from __future__ import annotations
import numpy as np

# 1D linear FE heat conduction, theta-method
# nodes: (n,1) or (n,), elems: (ne,2), k,rho,cp, dt,t_end
def solve_heat1d_linear_transient(env):
    x = np.asarray(env.get('nodes'), float).reshape(-1)
    elems = np.asarray(env.get('elems'), int)
    k = float(env.get('k')); rho = float(env.get('rho')); cp = float(env.get('cp'))
    dt = float(env.get('dt')); t_end = float(env.get('t_end')); theta = float(env.get('theta', 1.0))
    q = env.get('q', None)
    T0 = env.get('T0', 0.0)

    nn = x.shape[0]; ndof = nn
    K = np.zeros((ndof, ndof)); M = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    for e,(n1,n2) in enumerate(elems):
        x1,x2 = x[n1], x[n2]
        L = float(abs(x2-x1))
        if L <= 0:
            raise ValueError("Zero-length 1D element")
        Ke = (k/L)*np.array([[1,-1],[-1,1]], float)
        Me = (rho*cp*L/6.0)*np.array([[2,1],[1,2]], float)
        Fe = np.zeros(2)
        if q is not None:
            qv = float(q)
            Fe += qv * L/2.0 * np.array([1.0,1.0], float)
        dofs = [int(n1), int(n2)]
        for i,I in enumerate(dofs):
            F[I] += Fe[i]
            for j,J in enumerate(dofs):
                K[I,J] += Ke[i,j]
                M[I,J] += Me[i,j]

    # Dirichlet BC
    fixed = {}
    for ent in env.get('fixT', []):
        if len(ent) == 2:
            n = int(ent[0]); val = float(ent[1]); fixed[n] = val
    fixed_idx = np.array(sorted(fixed.keys()), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed_idx)

    # initial condition
    if np.ndim(T0) == 0:
        T = np.zeros(ndof) + float(T0)
    else:
        T = np.asarray(T0, float).reshape(ndof)
    for i in fixed_idx:
        T[i] = fixed[i]

    nstep = int(np.ceil(t_end/dt))
    times = np.linspace(0.0, nstep*dt, nstep+1)
    Thist = np.zeros((nstep+1, ndof))
    Thist[0,:] = T

    A = M + theta*dt*K

    for n in range(nstep):
        b = (M - (1.0-theta)*dt*K) @ T + dt*F
        if fixed_idx.size > 0:
            b_eff = b[free] - A[np.ix_(free, fixed_idx)] @ np.array([fixed[i] for i in fixed_idx], float)
            Tfree = np.linalg.solve(A[np.ix_(free, free)], b_eff)
            T[free] = Tfree
            T[fixed_idx] = np.array([fixed[i] for i in fixed_idx], float)
        else:
            T = np.linalg.solve(A, b)
        Thist[n+1,:] = T

    return {'T_hist': Thist, 'times': times, 'nodes': x}
