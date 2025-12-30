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

from .mesh_rect import rect_mesh
import numpy as np

def solve_elas2d_plane_stress_thermal(env):
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float); elems = np.array(env["elems"], dtype=int)
    else:
        Lx=float(env.get("Lx",1.0)); Ly=float(env.get("Ly",1.0))
        nx=int(env.get("nx",12)); ny=int(env.get("ny",12))
        nodes, elems = rect_mesh(Lx,Ly,nx,ny)
    E=float(env["E"]); nu=float(env["nu"]); t=float(env.get("t",1.0))
    alpha=float(env.get("alpha", 0.0))
    D = (E/(1-nu**2))*np.array([[1,  nu, 0],
                                [nu,  1, 0],
                                [0,   0, (1-nu)/2]], dtype=float)
    nn=nodes.shape[0]
    if "dT" in env:
        dT = env["dT"]
        if isinstance(dT, (list, tuple, np.ndarray)) and len(dT)==nn:
            dT_n = np.array(dT, dtype=float)
        else:
            dT_n = np.full(nn, float(dT), dtype=float)
    else:
        dT_n = np.zeros(nn)

    ndof=2*nn; K=np.zeros((ndof,ndof)); F=np.zeros(ndof)
    qp,qw=_gauss2()
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
            Ke += (B.T@D@B) * t * detJ * wt
            dTe = N @ dT_n[conn]
            eps_th = alpha * dTe * np.array([1.0,1.0,0.0])
            Fe += - (B.T @ (D @ eps_th)) * t * detJ * wt
        dofs=[]
        for n in conn: dofs.extend([2*n,2*n+1])
        for i in range(8):
            for j in range(8):
                K[dofs[i],dofs[j]]+=Ke[i,j]
        for i in range(8): F[dofs[i]]+=Fe[i]

    fixed=[]
    if "supports" in env:
        for n,fx,fy in np.array(env["supports"], dtype=int):
            if fx: fixed.append(2*int(n))
            if fy: fixed.append(2*int(n)+1)
    xs=nodes[:,0]; ys=nodes[:,1]; tol=1e-12
    minx,maxx,miny,maxy = xs.min(),xs.max(),ys.min(),ys.max()
    for key,mask in [
        ("clamp_left",   np.abs(xs-minx)<tol),
        ("clamp_right",  np.abs(xs-maxx)<tol),
        ("clamp_bottom", np.abs(ys-miny)<tol),
        ("clamp_top",    np.abs(ys-maxy)<tol),
    ]:
        if int(env.get(key,0))==1:
            ids=np.where(mask)[0]
            for i in ids: fixed.extend([2*i,2*i+1])
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    U=np.zeros(ndof); U[free]=np.linalg.solve(K[free][:,free], F[free])

    stress=[]
    for conn in elems:
        xe=nodes[conn,:]
        xi=eta=0.0
        N,dN_dxi=_shape_Q4(xi,eta); J=dN_dxi@xe; invJ=np.linalg.inv(J); dN_dx=invJ@dN_dxi
        B=np.zeros((3,8))
        for i in range(4):
            B[0,2*i]   = dN_dx[0,i]
            B[1,2*i+1] = dN_dx[1,i]
            B[2,2*i]   = dN_dx[1,i]
            B[2,2*i+1] = dN_dx[0,i]
        dofs=[]; 
        for n in conn: dofs.extend([2*n,2*n+1])
        eps=B@U[dofs]; dTe = N @ dT_n[conn]
        sig=D@(eps - alpha*dTe*np.array([1.0,1.0,0.0]))
        stress.append(sig)
    return {"nodes":nodes,"elems":elems,"U":U,"stress":np.array(stress)}
