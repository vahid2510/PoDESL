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

def _D_plane_stress(E, nu):
    c = E/(1-nu**2)
    return c*np.array([[1, nu, 0],
                       [nu, 1, 0],
                       [0,  0, (1-nu)/2]])

def solve_elas2d_planestress(env):
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float)
        elems = np.array(env["elems"], dtype=int)
    else:
        Lx = float(env.get("Lx", env.get("L", 1.0)))
        Ly = float(env.get("Ly", 1.0))
        nx = int(env.get("nx", 10)); ny = int(env.get("ny", 10))
        nodes, elems = rect_mesh(Lx, Ly, nx, ny)

    t = float(env.get("t", 1.0))
    E = float(env["E"]); nu = float(env["nu"])
    D = _D_plane_stress(E, nu)

    nn = nodes.shape[0]; ndof = nn*2
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    qp, qw = _gauss2()
    for e in range(elems.shape[0]):
        conn = elems[e]; xe = nodes[conn,:]
        Ke = np.zeros((8,8)); Fe = np.zeros(8)
        for (xi,eta), wt in zip(qp, qw):
            N, dN_dxi = _shape_Q4(xi, eta)
            J = dN_dxi @ xe; detJ = np.linalg.det(J); invJ = np.linalg.inv(J)
            dN_dx = invJ @ dN_dxi
            B = np.zeros((3,8))
            for i in range(4):
                B[0, 2*i+0] = dN_dx[0,i]
                B[1, 2*i+1] = dN_dx[1,i]
                B[2, 2*i+0] = dN_dx[1,i]
                B[2, 2*i+1] = dN_dx[0,i]
            Ke += (B.T @ D @ B) * t * detJ * wt
            bx = float(env.get("bx", 0.0)); by = float(env.get("by", 0.0))
            fe_xy = np.zeros(8)
            for i in range(4):
                fe_xy[2*i+0] += N[i] * bx * t * detJ * wt
                fe_xy[2*i+1] += N[i] * by * t * detJ * wt
            Fe += fe_xy
        dofs = []
        for n in conn:
            dofs.extend([2*n+0, 2*n+1])
        for i in range(8):
            I = dofs[i]; F[I] += Fe[i]
            for j in range(8):
                J = dofs[j]
                K[I,J] += Ke[i,j]

    if "loads" in env:
        loads = np.array(env["loads"], dtype=float)
        for ld in loads:
            n, Fx, Fy = int(ld[0]), float(ld[1]), float(ld[2])
            F[2*n+0] += Fx; F[2*n+1] += Fy

    fixed = []
    if "supports" in env:
        supports = np.array(env["supports"], dtype=float)
        for s in supports:
            n, fx, fy = int(s[0]), int(s[1]), int(s[2])
            if fx: fixed.append(2*n+0)
            if fy: fixed.append(2*n+1)

    xs = nodes[:,0]; ys = nodes[:,1]
    tol = 1e-12
    def clamp_edge(flag, mask):
        v = env.get(flag, None)
        if v is None: return
        if int(v) == 1:
            ids = np.where(mask)[0]
            for i in ids: 
                fixed.append(2*i+0); fixed.append(2*i+1)

    clamp_edge("clamp_left",   np.abs(xs - xs.min()) < tol)
    clamp_edge("clamp_right",  np.abs(xs - xs.max()) < tol)
    clamp_edge("clamp_bottom", np.abs(ys - ys.min()) < tol)
    clamp_edge("clamp_top",    np.abs(ys - ys.max()) < tol)

    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed)
    if free.size == 0:
        raise ValueError("No free DOFs. Check supports/edge clamps.")
    U = np.zeros(ndof)
    U[free] = np.linalg.solve(K[free][:, free], F[free])

    # stress at centroid
    stress = []
    for e in range(elems.shape[0]):
        conn = elems[e]; xe = nodes[conn,:]
        N, dN_dxi = _shape_Q4(0.0, 0.0)
        J = dN_dxi @ xe; invJ = np.linalg.inv(J)
        dN_dx = invJ @ dN_dxi
        B = np.zeros((3,8)); dofs = []
        for i in range(4):
            B[0, 2*i+0] = dN_dx[0,i]
            B[1, 2*i+1] = dN_dx[1,i]
            B[2, 2*i+0] = dN_dx[1,i]
            B[2, 2*i+1] = dN_dx[0,i]
            dofs.extend([2*conn[i]+0, 2*conn[i]+1])
        ue = U[dofs]
        sig = D @ (B @ ue)
        stress.append(sig)
    stress = np.array(stress)

    return {"nodes": nodes, "elems": elems, "U": U, "stress": stress}
