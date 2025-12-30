import numpy as np
from .mesh_rect import rect_mesh

def _gauss2():
    a = 1/np.sqrt(3); pts = [(-a,-a),(a,-a),(a,a),(-a,a)]; w=[1,1,1,1]; return pts,w
def _shape_Q4(xi,eta):
    N = 0.25*np.array([(1-xi)*(1-eta),
                       (1+xi)*(1-eta),
                       (1+xi)*(1+eta),
                       (1-xi)*(1+eta)], dtype=float)
    dN_dxi = 0.25*np.array([[-(1-eta),  (1-eta), (1+eta), -(1+eta)],
                            [-(1-xi),  -(1+xi), (1+xi),   (1-xi)]], dtype=float)
    return N, dN_dxi

def solve_heat2d_transient_conv(env):
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float); elems = np.array(env["elems"], dtype=int)
    else:
        Lx=float(env.get("Lx",1.0)); Ly=float(env.get("Ly",1.0))
        nx=int(env.get("nx",20)); ny=int(env.get("ny",20))
        nodes, elems = rect_mesh(Lx,Ly,nx,ny)
    kx=float(env.get("kx", 30.0)); ky=float(env.get("ky", 30.0))
    rho=float(env.get("rho", 7850.0)); cp=float(env.get("cp", 500.0))
    dt=float(env.get("dt", 0.01)); t_end=float(env.get("t_end", 1.0)); theta=float(env.get("theta", 1.0))
    T0 = env.get("T0", 20.0)
    qdot=float(env.get("qdot", 0.0))
    hL=float(env.get("h_left", 0.0));  TinfL=env.get("Tinf_left", 0.0)
    hR=float(env.get("h_right",0.0));  TinfR=env.get("Tinf_right",0.0)
    hB=float(env.get("h_bottom",0.0)); TinfB=env.get("Tinf_bottom",0.0)
    hT=float(env.get("h_top",0.0));    TinfT=env.get("Tinf_top",0.0)

    nn=nodes.shape[0]; ndof=nn
    K=np.zeros((ndof,ndof)); C=np.zeros((ndof,ndof)); F=np.zeros(ndof)
    qp,qw=_gauss2()
    for conn in elems:
        xe=nodes[conn,:]
        Ke=np.zeros((4,4)); Ce=np.zeros((4,4)); Fe=np.zeros(4)
        for (xi,eta),wt in zip(qp,qw):
            N,dN_dxi=_shape_Q4(xi,eta)
            J=dN_dxi@xe; detJ=np.linalg.det(J); invJ=np.linalg.inv(J)
            dN_dx=invJ@dN_dxi
            k_mat = np.diag([kx, ky])
            B = dN_dx
            Ke += (B.T@k_mat@B) * detJ * wt
            Ce += (np.outer(N,N) * rho*cp) * detJ * wt
            Fe += (N * qdot) * detJ * wt
        dofs=conn
        for i in range(4):
            for j in range(4):
                K[dofs[i],dofs[j]]+=Ke[i,j]
                C[dofs[i],dofs[j]]+=Ce[i,j]
        for i in range(4): F[dofs[i]]+=Fe[i]

    xs=nodes[:,0]; ys=nodes[:,1]; tol=1e-12
    minx,maxx,miny,maxy = xs.min(),xs.max(),ys.min(),ys.max()
    def add_conv(ids, h, Tinf):
        if h==0: return
        ids = ids[np.argsort(ids)]
        ids = ids[np.argsort(ys[ids] if np.allclose(xs[ids], xs[ids][0]) else xs[ids])]
        for a,b in zip(ids[:-1], ids[1:]):
            xa,ya = nodes[a]; xb,yb = nodes[b]
            Le = np.hypot(xb-xa,yb-ya)
            Ke = h*Le/6.0 * np.array([[2,1],[1,2]], dtype=float)
            Fe = h*Le/2.0 * Tinf * np.array([1,1], dtype=float)
            K[a,a]+=Ke[0,0]; K[a,b]+=Ke[0,1]; K[b,a]+=Ke[1,0]; K[b,b]+=Ke[1,1]
            F[a]+=Fe[0]; F[b]+=Fe[1]
    left_ids   = np.where(np.abs(xs-minx)<tol)[0]
    right_ids  = np.where(np.abs(xs-maxx)<tol)[0]
    bottom_ids = np.where(np.abs(ys-miny)<tol)[0]
    top_ids    = np.where(np.abs(ys-maxy)<tol)[0]
    add_conv(left_ids, hL, TinfL)
    add_conv(right_ids, hR, TinfR)
    add_conv(bottom_ids, hB, TinfB)
    add_conv(top_ids, hT, TinfT)

    if isinstance(T0, (list, tuple, np.ndarray)) and len(T0)==nn:
        T = np.array(T0, dtype=float)
    else:
        T = np.full(nn, float(T0), dtype=float)

    nsteps=int(np.ceil(t_end/dt))
    times=np.linspace(0,nsteps*dt,nsteps+1)
    A_theta = C + theta*dt*K
    bfac1 = C - (1-theta)*dt*K
    Thist=[T.copy()]
    for step in range(1, nsteps+1):
        rhs = bfac1@T + dt*F
        T = np.linalg.solve(A_theta, rhs)
        Thist.append(T.copy())
    return {"times":times, "T_hist":np.array(Thist)}
