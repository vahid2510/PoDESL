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

def _assemble_domain(nodes, elems, k=1.0, t=1.0, q=None):
    nn = nodes.shape[0]
    K = np.zeros((nn, nn)); F = np.zeros(nn)
    qp, qw = _gauss2()
    for e in range(elems.shape[0]):
        conn = elems[e]; xe = nodes[conn,:]
        Ke = np.zeros((4,4)); Fe = np.zeros(4)
        for (xi,eta), wt in zip(qp, qw):
            N, dN_dxi = _shape_Q4(xi, eta)
            J = dN_dxi @ xe; detJ = np.linalg.det(J); invJ = np.linalg.inv(J)
            dN_dx = invJ @ dN_dxi
            Ke += (k*t) * (dN_dx.T @ dN_dx) * detJ * wt
            if q is not None:
                Fe += N * (q * t) * detJ * wt
        for i in range(4):
            I = conn[i]; F[I] += Fe[i]
            for j in range(4):
                J = conn[j]; K[I,J] += Ke[i,j]
    return K, F

def _edge_gauss2():
    a = 1/np.sqrt(3)
    return [(-a,1),(a,1)], [1,1]  # param for edge with mapping, we will override

def _apply_edge_bc(nodes, elems, K, F, k=1.0, t=1.0, flags=None):
    # flags: dict of possible keys:
    # - Dirichlet: T_left, T_right, T_top, T_bottom
    # - Neumann:   qn_left, qn_right, qn_top, qn_bottom  (W/m^2) positive outward
    # - Robin:     h_left, Tinf_left  (and right/top/bottom)
    xs = nodes[:,0]; ys = nodes[:,1]; tol = 1e-12
    minx = xs.min(); maxx = xs.max(); miny = ys.min(); maxy = ys.max()
    N1D = [lambda s: 0.5*(1-s), lambda s: 0.5*(1+s)]
    s_pts = [-1/np.sqrt(3), 1/np.sqrt(3)]
    w_pts = [1,1]

    # collect edges as sequences of node ids per side (rectangular structured assumed)
    left_ids   = np.where(np.abs(xs - minx) < tol)[0]
    right_ids  = np.where(np.abs(xs - maxx) < tol)[0]
    bottom_ids = np.where(np.abs(ys - miny) < tol)[0]
    top_ids    = np.where(np.abs(ys - maxy) < tol)[0]

    def sort_by_y(ids): return ids[np.argsort(ys[ids])]
    def sort_by_x(ids): return ids[np.argsort(xs[ids])]

    left_ids, right_ids   = sort_by_y(left_ids), sort_by_y(right_ids)
    bottom_ids, top_ids   = sort_by_x(bottom_ids), sort_by_x(top_ids)

    # helper to integrate along polyline edges (between consecutive nodes)
    def edge_contrib(ids, h_coef=None, Tinf=None, qn=None):
        for a,b in zip(ids[:-1], ids[1:]):
            xa, ya = xs[a], ys[a]
            xb, yb = xs[b], ys[b]
            L = np.hypot(xb-xa, yb-ya)
            # 2-node line shape functions
            for sp, wp in zip(s_pts, w_pts):
                N = np.array([N1D[0](sp), N1D[1](sp)])
                J = L/2.0
                # Map to global DOF indices
                conn = [a,b]
                if qn is not None:
                    Fe = N * (qn * t) * J * wp
                    F[conn[0]] += Fe[0]; F[conn[1]] += Fe[1]
                if h_coef is not None and Tinf is not None:
                    Ke = (h_coef * t) * (np.outer(N,N)) * J * wp
                    Fe = (h_coef * Tinf * t) * N * J * wp
                    for i in range(2):
                        F[conn[i]] += Fe[i]
                        for j in range(2):
                            K[conn[i], conn[j]] += Ke[i,j]

    # Dirichlet sets
    fixed = []
    fixed_vals = {}

    # Dirichlet
    if "T_left" in flags:
        for i in left_ids: fixed.append(i); fixed_vals[i] = float(flags["T_left"])
    if "T_right" in flags:
        for i in right_ids: fixed.append(i); fixed_vals[i] = float(flags["T_right"])
    if "T_bottom" in flags:
        for i in bottom_ids: fixed.append(i); fixed_vals[i] = float(flags["T_bottom"])
    if "T_top" in flags:
        for i in top_ids: fixed.append(i); fixed_vals[i] = float(flags["T_top"])

    # Neumann flux
    if "qn_left" in flags:   edge_contrib(left_ids,  qn=float(flags["qn_left"]))
    if "qn_right" in flags:  edge_contrib(right_ids, qn=float(flags["qn_right"]))
    if "qn_bottom" in flags: edge_contrib(bottom_ids, qn=float(flags["qn_bottom"]))
    if "qn_top" in flags:    edge_contrib(top_ids, qn=float(flags["qn_top"]))

    # Robin (convection)
    for side in ["left","right","bottom","top"]:
        hkey = f"h_{side}"; Tkey = f"Tinf_{side}"
        if hkey in flags and Tkey in flags:
            ids = {"left":left_ids,"right":right_ids,"bottom":bottom_ids,"top":top_ids}[side]
            edge_contrib(ids, h_coef=float(flags[hkey]), Tinf=float(flags[Tkey]))

    return fixed, fixed_vals

def solve_heat2d_steady_bc(env):
    # mesh
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float)
        elems = np.array(env["elems"], dtype=int)
    else:
        Lx = float(env.get("Lx", 1.0)); Ly = float(env.get("Ly", 1.0))
        nx = int(env.get("nx", 20)); ny = int(env.get("ny", 20))
        nodes, elems = rect_mesh(Lx, Ly, nx, ny)
    k = float(env.get("k", 1.0)); t = float(env.get("t", 1.0))
    q = env.get("q", None); qval = float(q) if q is not None else None

    K, F = _assemble_domain(nodes, elems, k=k, t=t, q=qval)
    flags = {k:v for k,v in env.items() if isinstance(k,str)}
    fixed, fixed_vals = _apply_edge_bc(nodes, elems, K, F, k=k, t=t, flags=flags)

    # apply dirichlet
    A = K.copy(); b = F.copy()
    for i in fixed:
        A[i,:] = 0.0; A[:,i] = 0.0; A[i,i] = 1.0
        b[i] = fixed_vals[i]
    T = np.linalg.solve(A, b)
    return {"nodes": nodes, "elems": elems, "T": T}
