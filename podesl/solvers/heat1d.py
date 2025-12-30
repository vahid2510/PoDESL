import numpy as np

def solve_heat1d_steady(env):
    L = float(env.get("L"))
    k = float(env.get("k"))
    A = float(env.get("A", 1.0))
    qdot = float(env.get("qdot", 0.0))
    nel = int(env.get("nel", 10))
    left = str(env.get("left", "Dirichlet")).lower()
    right = str(env.get("right", "Dirichlet")).lower()

    nn = nel + 1
    h = L/nel
    K = np.zeros((nn, nn))
    f = np.zeros(nn)

    # element stiffness and load for linear 1D conduction
    ke = (k*A/h) * np.array([[1,-1],[-1,1]], dtype=float)
    feq = (qdot*A*h/2.0) * np.array([1,1], dtype=float)

    for e in range(nel):
        dofs = [e, e+1]
        for i in range(2):
            f[dofs[i]] += feq[i]
            for j in range(2):
                K[dofs[i], dofs[j]] += ke[i, j]

    # boundary at left
    if left == "dirichlet" or left == "d":
        T_left = float(env.get("T_left"))
        # enforce by big penalty or elimination; do elimination:
        K[0,:] = 0.0; K[:,0] = 0.0; K[0,0] = 1.0
        f[0] = T_left
    elif left == "robin" or left == "h":
        hL = float(env.get("h_left")); TinfL = float(env.get("Tinf_left"))
        K[0,0] += hL*A
        f[0] += hL*A*TinfL
    else:
        pass

    # boundary at right
    if right == "dirichlet" or right == "d":
        T_right = float(env.get("T_right"))
        K[-1,:] = 0.0; K[:,-1] = 0.0; K[-1,-1] = 1.0
        f[-1] = T_right
    elif right == "robin" or right == "h":
        hR = float(env.get("h_right")); TinfR = float(env.get("Tinf_right"))
        K[-1,-1] += hR*A
        f[-1] += hR*A*TinfR
    else:
        pass

    T = np.linalg.solve(K, f)
    x = np.linspace(0.0, L, nn)
    return {"T": T, "x": x}
