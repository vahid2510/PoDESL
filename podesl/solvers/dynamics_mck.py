import numpy as np

def solve_mck_newmark(env):
    M = np.array(env["M"], dtype=float)
    C = np.array(env.get("C", np.zeros_like(M)), dtype=float)
    K = np.array(env["K"], dtype=float)
    dt = float(env["dt"]); t_end = float(env["t_end"])
    u0 = np.array(env.get("u0", np.zeros(M.shape[0])), dtype=float)
    v0 = np.array(env.get("v0", np.zeros(M.shape[0])), dtype=float)

    n = M.shape[0]
    nt = int(round(t_end/dt)) + 1
    t = np.linspace(0.0, t_end, nt)

    f_env = env.get("f", np.zeros((nt, n)))
    if isinstance(f_env, (list, tuple, np.ndarray)):
        f = np.array(f_env, dtype=float)
        if f.ndim == 1:
            f = np.tile(f, (nt, 1))
        elif f.shape[0] != nt:
            if f.shape[0] < nt:
                pad = np.tile(f[-1], (nt - f.shape[0], 1))
                f = np.vstack([f, pad])
            else:
                f = f[:nt, :]
    else:
        f = np.zeros((nt, n))

    beta = 1/4; gamma = 1/2
    A0 = 1.0/(beta*dt*dt); A1 = gamma/(beta*dt)
    Keff = K + A1*C + A0*M

    chol = None
    if np.allclose(Keff, Keff.T, atol=1e-12):
        try:
            chol = np.linalg.cholesky(Keff)
        except Exception:
            chol = None

    def solve_keff(rhs):
        if chol is not None:
            y = np.linalg.solve(chol, rhs)
            return np.linalg.solve(chol.T, y)
        return np.linalg.solve(Keff, rhs)

    U = np.zeros((nt, n)); V = np.zeros_like(U); A = np.zeros_like(U)
    U[0] = u0; V[0] = v0
    A[0] = np.linalg.solve(M, f[0] - C@V[0] - K@U[0])

    for k in range(nt-1):
        u_pred = U[k] + dt*V[k] + dt*dt*(0.5 - beta)*A[k]
        v_pred = V[k] + dt*(1 - gamma)*A[k]
        rhs = f[k+1] + M.dot(A0*u_pred) + C.dot(A1*u_pred + v_pred)
        U[k+1] = solve_keff(rhs)
        A[k+1] = A0*(U[k+1] - u_pred)
        V[k+1] = v_pred + gamma*dt*A[k+1]

    return {"t": t, "U": U, "V": V, "A": A}
