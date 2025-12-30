import numpy as np
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
def solve_thermomech2d_ps_coupled(env):
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float); elems = np.array(env["elems"], dtype=int)
    else:
        Lx=float(env.get("Lx",1.0)); Ly=float(env.get("Ly",1.0))
        nx=int(env.get("nx",10)); ny=int(env.get("ny",10))
        xs = np.linspace(0.0, Lx, nx+1)
        ys = np.linspace(0.0, Ly, ny+1)
        nodes = np.array([[x,y] for j in range(ny+1) for x in xs for y in [ys[j]]], dtype=float)
        def nid(i,j): return j*(nx+1) + i
        elems = []
        for j in range(ny):
            for i in range(nx):
                n00 = nid(i, j); n10 = nid(i+1, j); n11 = nid(i+1, j+1); n01 = nid(i, j+1)
                elems.append([n00, n10, n11, n01])
        elems = np.array(elems, dtype=int)
    E=float(env["E"]); nu=float(env["nu"]); t=float(env.get("t",1.0))
    alpha=float(env["alpha"])
    dT = env.get("dT", None); T = env.get("T", None); Tref = float(env.get("Tref", 0.0))
    D = (E/(1-nu**2))*np.array([[1, nu, 0],[nu, 1, 0],[0,  0, (1-nu)/2]], dtype=float)
    nn=nodes.shape[0]; ndof=2*nn
    K=np.zeros((ndof,ndof)); F=np.zeros(ndof)
    qp,qw=_gauss2()
    if T is not None:
        Tn = np.array(T, dtype=float).reshape(-1)
        if Tn.size!=nn: raise ValueError("Length of T must equal nnodes")
        dTn = Tn - Tref
    elif dT is not None:
        dTn = np.full(nn, float(dT), dtype=float)
    else:
        dTn = np.zeros(nn)
    for conn in elems:
        xe=nodes[conn,:]
        Ke=np.zeros((8,8)); Fe=np.zeros(8)
        for (xi,eta),wt in zip(qp,qw):
            N,dN_dxi=_shape_Q4(xi,eta)
            J=dN_dxi@xe; detJ=np.linalg.det(J); invJ=np.linalg.inv(J)
            dN_dx=invJ@dN_dxi
            B=np.zeros((3,8))
            for i in range(4):
                B[0,2*i]   = dN_dx[0,i]
                B[1,2*i+1] = dN_dx[1,i]
                B[2,2*i]   = dN_dx[1,i]
                B[2,2*i+1] = dN_dx[0,i]
            Ke += (B.T@D@B) * detJ*wt * t
            dT_gp = float(np.dot(N, dTn[conn]))
            eps_th = alpha*dT_gp*np.array([1,1,0], dtype=float)
            Fe += (B.T@D@eps_th) * detJ*wt * t
        dofs=[]
        for n in conn: dofs.extend([2*n,2*n+1])
        for i in range(8):
            for j in range(8):
                K[dofs[i],dofs[j]]+=Ke[i,j]
        for i in range(8): F[dofs[i]]+=Fe[i]
    fixed=[]
    for n,fx,fy in env.get("supports", []):
        n=int(n)
        if fx: fixed.append(2*n)
        if fy: fixed.append(2*n+1)
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    U=np.zeros(ndof); U[free]=np.linalg.solve(K[free][:,free], F[free])
    return {"U":U}
