import numpy as np
from .mesh_rect import rect_mesh

def _gauss2():
    a = 1/np.sqrt(3)
    pts = [(-a,-a),(a,-a),(a,a),(-a,a)]
    w = [1,1,1,1]
    return pts, w

def _shape_Q4(xi, eta):
    N = 0.25*np.array([(1-xi)*(1-eta),
                       (1+xi)*(1-eta),
                       (1+xi)*(1+eta),
                       (1-xi)*(1+eta)], dtype=float)
    dN_dxi = 0.25*np.array([[-(1-eta),  (1-eta), (1+eta), -(1+eta)],
                            [-(1-xi),  -(1+xi), (1+xi),   (1-xi)]], dtype=float)
    return N, dN_dxi

def _assemble(nodes, elems, k=1.0, t=1.0, q=None):
    nn = nodes.shape[0]
    K = np.zeros((nn, nn)); F = np.zeros(nn)
    qp, qw = _gauss2()
    for e in range(elems.shape[0]):
        conn = elems[e]; xe = nodes[conn,:]
        Ke = np.zeros((4,4)); Fe = np.zeros(4)
        for (xi,eta), wt in zip(qp, qw):
            N, dN_dxi = _shape_Q4(xi, eta)
            J = dN_dxi @ xe
            detJ = np.linalg.det(J); invJ = np.linalg.inv(J)
            dN_dx = invJ @ dN_dxi
            B = dN_dx
            Ke += (k*t) * (B.T @ B) * detJ * wt
            if q is not None:
                Fe += N * (q * t) * detJ * wt
        for i in range(4):
            I = conn[i]
            F[I] += Fe[i]
            for j in range(4):
                J = conn[j]
                K[I,J] += Ke[i,j]
    return K, F

def solve_heat2d_steady(env):
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float)
        elems = np.array(env["elems"], dtype=int)
    else:
        Lx = float(env.get("Lx", env.get("L", 1.0)))
        Ly = float(env.get("Ly", 1.0))
        nx = int(env.get("nx", 10)); ny = int(env.get("ny", 10))
        nodes, elems = rect_mesh(Lx, Ly, nx, ny)
    k = float(env.get("k", 1.0)); t = float(env.get("t", 1.0))
    q = env.get("q", None)
    qval = float(q) if q is not None else None

    K, F = _assemble(nodes, elems, k=k, t=t, q=qval)

    xs = nodes[:,0]; ys = nodes[:,1]
    tol = 1e-12
    fixed = []; fixed_vals = {}

    if "T_left" in env:
        ids = np.where(np.abs(xs - xs.min()) < tol)[0]
        for i in ids: fixed.append(i); fixed_vals[i] = float(env["T_left"])
    if "T_right" in env:
        ids = np.where(np.abs(xs - xs.max()) < tol)[0]
        for i in ids: fixed.append(i); fixed_vals[i] = float(env["T_right"])
    if "T_bottom" in env:
        ids = np.where(np.abs(ys - ys.min()) < tol)[0]
        for i in ids: fixed.append(i); fixed_vals[i] = float(env["T_bottom"])
    if "T_top" in env:
        ids = np.where(np.abs(ys - ys.max()) < tol)[0]
        for i in ids: fixed.append(i); fixed_vals[i] = float(env["T_top"])

    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(nodes.shape[0]), fixed)

    T = np.zeros(nodes.shape[0])
    A = K.copy(); b = F.copy()
    for i in fixed:
        A[i,:] = 0.0; A[:,i] = 0.0; A[i,i] = 1.0
        b[i] = fixed_vals[i]
    if free.size>0:
        T = np.linalg.solve(A, b)
    else:
        T = b

    return {"nodes": nodes, "elems": elems, "T": T}
