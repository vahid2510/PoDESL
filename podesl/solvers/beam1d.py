import numpy as np

def _k_e(EI, h):
    a = EI / h**3
    ke = a * np.array([
        [ 12,   6*h,  -12,   6*h],
        [ 6*h, 4*h*h, -6*h, 2*h*h],
        [-12,  -6*h,  12,  -6*h],
        [ 6*h, 2*h*h, -6*h, 4*h*h],
    ], dtype=float)
    return ke

def _uniform_load_vector(q, h):
    return q * h / 12.0 * np.array([6, 3*h, 6, -3*h], dtype=float)

def solve_beam_static(env):
    L = float(env.get("L"))
    E = float(env.get("E"))
    I = float(env.get("I"))
    nel = int(env.get("nel", 10))
    q = float(env.get("q", 0.0))
    point_loads = env.get("point_loads", []) or []
    left = str(env.get("left", "clamped")).lower()
    right = str(env.get("right", "clamped")).lower()

    nn = nel + 1
    ndof = nn*2
    h = L/nel
    K = np.zeros((ndof, ndof))
    f = np.zeros(ndof)

    EI = E*I
    for e in range(nel):
        ke = _k_e(EI, h)
        fe = _uniform_load_vector(q, h)
        i = e
        dofs = [2*i, 2*i+1, 2*i+2, 2*i+3]
        for a in range(4):
            f[dofs[a]] += fe[a]
            for b in range(4):
                K[dofs[a], dofs[b]] += ke[a, b]

    for P, xpos in point_loads:
        idx = int(round(xpos / h))
        idx = max(0, min(nn-1, idx))
        f[2*idx] += float(P)

    bc_dofs = []
    if left == "clamped": bc_dofs += [0,1]
    elif left == "pinned": bc_dofs += [0]
    if right == "clamped": bc_dofs += [ndof-2, ndof-1]
    elif right == "pinned": bc_dofs += [ndof-2]

    free = np.setdiff1d(np.arange(ndof), np.array(bc_dofs, dtype=int))
    if free.size == 0: raise ValueError("No free DOFs. Check boundary conditions.")

    Kff = K[free][:, free]
    ff = f[free]
    u = np.zeros(ndof)
    u_free = np.linalg.solve(Kff, ff)
    u[free] = u_free

    reactions = K @ u - f
    R = reactions[bc_dofs]

    w = u[0::2]; theta = u[1::2]
    M = np.zeros(nn)
    for i in range(1, nn-1):
        M[i] = EI * (w[i-1] - 2*w[i] + w[i+1]) / (h**2)

    x = np.linspace(0, L, nn)
    return {"w": w, "theta": theta, "M": M, "x": x, "u": u, "reactions": R, "bc_dofs": bc_dofs}
