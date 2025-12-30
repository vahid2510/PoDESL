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

def solve_elas2d_ps_edgeload(env):
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float); elems = np.array(env["elems"], dtype=int)
    else:
        Lx = float(env.get("Lx", 1.0)); Ly = float(env.get("Ly", 1.0))
        nx = int(env.get("nx", 12)); ny = int(env.get("ny", 6))
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
        # assemble
        dofs = []
        for n in conn: dofs.extend([2*n+0, 2*n+1])
        for i in range(8):
            for j in range(8):
                K[dofs[i], dofs[j]] += Ke[i,j]

    # nodal loads
    if "loads" in env:
        loads = np.array(env["loads"], dtype=float)
        for ld in loads:
            n, Fx, Fy = int(ld[0]), float(ld[1]), float(ld[2])
            F[2*n+0] += Fx; F[2*n+1] += Fy

    # distributed edge tractions (rectangle edges only, constant)
    # tx_left, ty_left, ... act along edge length, contribute to nodal F
    xs = nodes[:,0]; ys = nodes[:,1]; tol = 1e-12
    minx, maxx, miny, maxy = xs.min(), xs.max(), ys.min(), ys.max()
    def sort_by_y(ids): return ids[np.argsort(ys[ids])]
    def sort_by_x(ids): return ids[np.argsort(xs[ids])]
    left_ids = sort_by_y(np.where(np.abs(xs-minx)<tol)[0])
    right_ids= sort_by_y(np.where(np.abs(xs-maxx)<tol)[0])
    bottom_ids=sort_by_x(np.where(np.abs(ys-miny)<tol)[0])
    top_ids   =sort_by_x(np.where(np.abs(ys-maxy)<tol)[0])

    N1 = lambda s: 0.5*(1-s); N2 = lambda s: 0.5*(1+s)
    s_pts = [-1/np.sqrt(3), 1/np.sqrt(3)]; w_pts=[1,1]

    def edge_load(ids, tx, ty):
        if tx==0.0 and ty==0.0: return
        for a,b in zip(ids[:-1], ids[1:]):
            xa,ya = xs[a],ys[a]; xb,yb = xs[b],ys[b]
            L = np.hypot(xb-xa,yb-ya)
            for sp,wp in zip(s_pts,w_pts):
                N = np.array([N1(sp), N2(sp)]); J = L/2.0
                fe = N * t * J * wp
                F[2*a+0] += fe[0]*tx; F[2*a+1] += fe[0]*ty
                F[2*b+0] += fe[1]*tx; F[2*b+1] += fe[1]*ty

    txl=float(env.get("tx_left",0.0)); tyl=float(env.get("ty_left",0.0))
    txr=float(env.get("tx_right",0.0)); tyr=float(env.get("ty_right",0.0))
    txb=float(env.get("tx_bottom",0.0)); tyb=float(env.get("ty_bottom",0.0))
    txt=float(env.get("tx_top",0.0)); tyt=float(env.get("ty_top",0.0))
    edge_load(left_ids, txl, tyl); edge_load(right_ids, txr, tyr)
    edge_load(bottom_ids, txb, tyb); edge_load(top_ids, txt, tyt)

    # supports
    fixed = []
    if "supports" in env:
        supports = np.array(env["supports"], dtype=float)
        for s in supports:
            n, fx, fy = int(s[0]), int(s[1]), int(s[2])
            if fx: fixed.append(2*n+0)
            if fy: fixed.append(2*n+1)
    for flag, ids in {"clamp_left":left_ids,"clamp_right":right_ids,"clamp_bottom":bottom_ids,"clamp_top":top_ids}.items():
        v = env.get(flag, None)
        if v is not None and int(v)==1:
            for i in ids: fixed.extend([2*i+0,2*i+1])

    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports/clamps.")
    U = np.zeros(ndof); U[free] = np.linalg.solve(K[free][:,free], F[free])

    # simple element stress at centroid (plane stress)
    stress = []
    for e in range(elems.shape[0]):
        conn = elems[e]; xe = nodes[conn,:]
        N, dN_dxi = _shape_Q4(0.0,0.0)
        J = dN_dxi @ xe; invJ = np.linalg.inv(J)
        dN_dx = invJ @ dN_dxi
        B = np.zeros((3,8)); dofs = []
        for i in range(4):
            B[0,2*i+0]=dN_dx[0,i]; B[1,2*i+1]=dN_dx[1,i]
            B[2,2*i+0]=dN_dx[1,i]; B[2,2*i+1]=dN_dx[0,i]
            dofs.extend([2*conn[i]+0, 2*conn[i]+1])
        ue = U[dofs]
        sig = D @ (B @ ue)
        stress.append(sig)
    stress = np.array(stress)

    return {"nodes": nodes, "elems": elems, "U": U, "stress": stress}
