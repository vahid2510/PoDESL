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

def solve_shell2d_mindlin_static(env):
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float); elems = np.array(env["elems"], dtype=int)
    else:
        Lx=float(env.get("Lx",1.0)); Ly=float(env.get("Ly",1.0))
        nx=int(env.get("nx",12)); ny=int(env.get("ny",12))
        nodes, elems = rect_mesh(Lx,Ly,nx,ny)
    E=float(env["E"]); nu=float(env["nu"]); t=float(env["t"])
    kappa=float(env.get("kappa", 5.0/6.0))
    q=float(env.get("q", 0.0))

    D_b = E*t**3/(12*(1-nu**2)) * np.array([[1, nu, 0],
                                            [nu, 1, 0],
                                            [0,  0, (1-nu)/2]], dtype=float)
    D_s = kappa*E*t/(2*(1+nu)) * np.eye(2)

    nn = nodes.shape[0]; ndof = 3*nn
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)
    qp,qw=_gauss2()

    for conn in elems:
        xe=nodes[conn,:]
        Ke=np.zeros((12,12)); Fe=np.zeros(12)
        for (xi,eta),wt in zip(qp,qw):
            N,dN_dxi=_shape_Q4(xi,eta)
            J=dN_dxi@xe; detJ=np.linalg.det(J); invJ=np.linalg.inv(J)
            dN_dx=invJ@dN_dxi
            Bb=np.zeros((3,12))
            for i in range(4):
                Bb[0, 3*i+1] = dN_dx[0,i]
                Bb[1, 3*i+2] = dN_dx[1,i]
                Bb[2, 3*i+1] = dN_dx[1,i]
                Bb[2, 3*i+2] = dN_dx[0,i]
            Bs=np.zeros((2,12))
            for i in range(4):
                Bs[0, 3*i]   = dN_dx[0,i]
                Bs[0, 3*i+1] = N[i]
                Bs[1, 3*i]   = dN_dx[1,i]
                Bs[1, 3*i+2] = N[i]
            Ke += (Bb.T@D_b@Bb + Bs.T@D_s@Bs) * detJ * wt
            Ne = np.zeros(12)
            for i in range(4): Ne[3*i] = N[i]
            Fe += (-q)*Ne * detJ * wt
        dofs=[]
        for n in conn: dofs.extend([3*n,3*n+1,3*n+2])
        for i in range(12):
            for j in range(12):
                K[dofs[i],dofs[j]]+=Ke[i,j]
        for i in range(12): F[dofs[i]]+=Fe[i]

    fixed=[]
    if "supports" in env:
        for n,fw,frx,fry in np.array(env["supports"], dtype=int):
            n=int(n)
            if fw:  fixed.append(3*n)
            if frx: fixed.append(3*n+1)
            if fry: fixed.append(3*n+2)
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
            for i in ids: fixed.extend([3*i,3*i+1,3*i+2])
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    U=np.zeros(ndof); U[free]=np.linalg.solve(K[free][:,free], F[free])
    return {"U":U}
