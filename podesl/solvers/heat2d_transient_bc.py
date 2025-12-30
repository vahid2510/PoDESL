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

def _assemble_MK(nodes, elems, rho, c, k, t, q=None):
    nn = nodes.shape[0]
    M = np.zeros((nn, nn)); K = np.zeros((nn, nn)); F = np.zeros(nn)
    qp, qw = _gauss2()
    for e in range(elems.shape[0]):
        conn = elems[e]; xe = nodes[conn,:]
        Me = np.zeros((4,4)); Ke = np.zeros((4,4)); Fe = np.zeros(4)
        for (xi,eta), wt in zip(qp, qw):
            N, dN_dxi = _shape_Q4(xi, eta)
            J = dN_dxi @ xe; detJ = np.linalg.det(J); invJ = np.linalg.inv(J)
            dN_dx = invJ @ dN_dxi
            Ke += (k * t) * (dN_dx.T @ dN_dx) * detJ * wt
            Me += (rho * c * t) * (np.outer(N, N)) * detJ * wt
            if q is not None:
                Fe += N * (q * t) * detJ * wt
        for i in range(4):
            I = conn[i]
            F[I] += Fe[i]
            for j in range(4):
                J = conn[j]
                M[I,J] += Me[i,j]
                K[I,J] += Ke[i,j]
    return M, K, F

def _interp_series(series, t):
    # series: [[t0, v0], [t1, v1], ...]; linear interpolation; outside => clamp
    arr = np.array(series, dtype=float)
    ts, vs = arr[:,0], arr[:,1]
    if t <= ts[0]: return vs[0]
    if t >= ts[-1]: return vs[-1]
    i = np.searchsorted(ts, t) - 1
    t0, t1 = ts[i], ts[i+1]; v0, v1 = vs[i], vs[i+1]
    a = (t - t0)/(t1 - t0)
    return (1-a)*v0 + a*v1

def _apply_edge_terms(nodes, K, F, t_thick, flags, now_t):
    # Neumann (qn_*), Robin (h_*, Tinf_*) possibly time dependent via *_series
    xs = nodes[:,0]; ys = nodes[:,1]
    tol = 1e-12
    minx,maxx,miny,maxy = xs.min(),xs.max(),ys.min(),ys.max()
    left   = np.where(np.abs(xs-minx)<tol)[0]
    right  = np.where(np.abs(xs-maxx)<tol)[0]
    bottom = np.where(np.abs(ys-miny)<tol)[0]
    top    = np.where(np.abs(ys-maxy)<tol)[0]

    def sort_y(ids): return ids[np.argsort(ys[ids])]
    def sort_x(ids): return ids[np.argsort(xs[ids])]
    left, right = sort_y(left), sort_y(right)
    bottom, top = sort_x(bottom), sort_x(top)

    N1 = lambda s: 0.5*(1-s); N2 = lambda s: 0.5*(1+s)
    s_pts = [-1/np.sqrt(3), 1/np.sqrt(3)]; w_pts=[1,1]

    def edge_pair(ids, qn=None, h=None, Tinf=None):
        for a,b in zip(ids[:-1], ids[1:]):
            xa,ya = xs[a],ys[a]; xb,yb = xs[b],ys[b]
            L = np.hypot(xb-xa,yb-ya); J = L/2.0
            for sp,wp in zip(s_pts,w_pts):
                N = np.array([N1(sp), N2(sp)])
                conn = [a,b]
                if qn is not None:
                    Fe = N * (qn * t_thick) * J * wp
                    F[conn[0]] += Fe[0]; F[conn[1]] += Fe[1]
                if (h is not None) and (Tinf is not None):
                    Ke = (h * t_thick) * np.outer(N,N) * J * wp
                    Fe = (h * Tinf * t_thick) * N * J * wp
                    for i in range(2):
                        F[conn[i]] += Fe[i]
                        for j in range(2):
                            K[conn[i], conn[j]] += Ke[i,j]

    # helper get value or series
    def val(name):
        if name in flags: return float(flags[name])
        sname = name + "_series"
        if sname in flags: return float(_interp_series(flags[sname], now_t))
        return None

    # Neumann
    edge_pair(left,   qn=val("qn_left"))
    edge_pair(right,  qn=val("qn_right"))
    edge_pair(bottom, qn=val("qn_bottom"))
    edge_pair(top,    qn=val("qn_top"))
    # Robin
    for side, ids in [("left",left),("right",right),("bottom",bottom),("top",top)]:
        h = val(f"h_{side}"); Tinf = val(f"Tinf_{side}")
        if (h is not None) and (Tinf is not None):
            edge_pair(ids, h=h, Tinf=Tinf)

def solve_heat2d_transient_bc(env):
    # mesh
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float); elems = np.array(env["elems"], dtype=int)
    else:
        Lx = float(env.get("Lx", 1.0)); Ly = float(env.get("Ly", 1.0))
        nx = int(env.get("nx", 20)); ny = int(env.get("ny", 20))
        nodes, elems = rect_mesh(Lx, Ly, nx, ny)

    rho = float(env["rho"]); c = float(env["c"]); k = float(env["k"]); t = float(env.get("t", 1.0))
    q = env.get("q", None); qval = float(q) if q is not None else None
    dt = float(env["dt"]); t_end = float(env["t_end"]); theta = float(env.get("theta", 0.5))
    nt = int(np.floor(t_end/dt)) + 1
    times = np.linspace(0.0, dt*(nt-1), nt)

    M, K, Fg = _assemble_MK(nodes, elems, rho, c, k, t, q=qval)

    nn = nodes.shape[0]
    T = np.zeros(nn)
    if "T_init" in env:
        Ti = env["T_init"]
        if isinstance(Ti, (list, tuple, np.ndarray)) and len(Ti) == nn:
            T = np.array(Ti, dtype=float)
        else:
            T.fill(float(Ti))

    xs = nodes[:,0]; ys = nodes[:,1]; tol = 1e-12
    fixed = []; fixed_vals = {}
    def set_edge(name, mask):
        if name in env:
            val = float(env[name]); ids = np.where(mask)[0]
            for i in ids: fixed.append(i); fixed_vals[i] = val
    set_edge("T_left",   np.abs(xs - xs.min()) < tol)
    set_edge("T_right",  np.abs(xs - xs.max()) < tol)
    set_edge("T_bottom", np.abs(ys - ys.min()) < tol)
    set_edge("T_top",    np.abs(ys - ys.max()) < tol)
    fixed = np.array(sorted(set(fixed)), dtype=int)

    A_imp_base = M/dt + theta*K
    A_exp_base = M/dt - (1.0-theta)*K

    Thist = np.zeros((nt, nn)); Thist[0,:] = T

    for n in range(1, nt):
        now_t = times[n]
        # copy and add edge terms for this time
        A_imp = A_imp_base.copy(); rhs = (A_exp_base @ T + Fg).copy()
        _apply_edge_terms(nodes, A_imp, rhs, t, env, now_t)

        # Dirichlet each step (kept constant over time here)
        A = A_imp; b = rhs
        for i in fixed:
            A[i,:] = 0.0; A[:,i] = 0.0; A[i,i] = 1.0
            b[i] = fixed_vals[i]
        T = np.linalg.solve(A, b)
        Thist[n,:] = T

    return {"nodes": nodes, "elems": elems, "t": times, "T": Thist, "T_end": Thist[-1,:]}
