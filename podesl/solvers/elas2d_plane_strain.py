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

def solve_elas2d_plane_strain(env):
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float); elems = np.array(env["elems"], dtype=int)
    else:
        Lx=float(env.get("Lx",1.0)); Ly=float(env.get("Ly",1.0))
        nx=int(env.get("nx",12)); ny=int(env.get("ny",12))
        nodes, elems = rect_mesh(Lx,Ly,nx,ny)
    E=float(env["E"]); nu=float(env["nu"]); t=float(env.get("t",1.0))
    c = E/((1+nu)*(1-2*nu))
    D = c*np.array([[1-nu,   nu,        0],
                    [  nu, 1-nu,        0],
                    [   0,    0, (1-2*nu)/2]], dtype=float)

    nn=nodes.shape[0]; ndof=2*nn
    K=np.zeros((ndof,ndof)); F=np.zeros(ndof)
    qp,qw=_gauss2()

    bx=float(env.get("bx",0.0)); by=float(env.get("by",0.0))

    for e,conn in enumerate(elems):
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
            Ne = np.zeros(8)
            for i in range(4):
                Ne[2*i]   = N[i]
                Ne[2*i+1] = N[i]
            Fe += (Ne * np.array([bx,by]*4)) * t * detJ * wt
        dofs=[]
        for n in conn: dofs.extend([2*n,2*n+1])
        for i in range(8):
            for j in range(8):
                K[dofs[i],dofs[j]]+=Ke[i,j]
        for i in range(8): F[dofs[i]]+=Fe[i]

    xs=nodes[:,0]; ys=nodes[:,1]; tol=1e-12
    minx,maxx,miny,maxy = xs.min(),xs.max(),ys.min(),ys.max()
    left=np.where(np.abs(xs-minx)<tol)[0]
    right=np.where(np.abs(xs-maxx)<tol)[0]
    bottom=np.where(np.abs(ys-miny)<tol)[0]
    top=np.where(np.abs(ys-maxy)<tol)[0]
    def sort_y(ids): return ids[np.argsort(ys[ids])]
    def sort_x(ids): return ids[np.argsort(xs[ids])]
    left,right = sort_y(left), sort_y(right)
    bottom,top = sort_x(bottom), sort_x(top)

    def add_edge(ids, tx=0.0, ty=0.0):
        N1=lambda s:0.5*(1-s); N2=lambda s:0.5*(1+s)
        s_pts=[-1/np.sqrt(3), 1/np.sqrt(3)]; w_pts=[1,1]
        for a,b in zip(ids[:-1], ids[1:]):
            xa,ya=xs[a],ys[a]; xb,yb=xs[b],ys[b]
            L=np.hypot(xb-xa,yb-ya); J=L/2.0
            for sp,wp in zip(s_pts,w_pts):
                N=np.array([N1(sp),N2(sp)])
                fe = np.array([tx,ty, tx,ty]) * J * wp * t
                dofs=[2*a,2*a+1,2*b,2*b+1]
                F[dofs]+= fe
    add_edge(left,   float(env.get("tx_left",0.0)),   float(env.get("ty_left",0.0)))
    add_edge(right,  float(env.get("tx_right",0.0)),  float(env.get("ty_right",0.0)))
    add_edge(bottom, float(env.get("tx_bottom",0.0)), float(env.get("ty_bottom",0.0)))
    add_edge(top,    float(env.get("tx_top",0.0)),    float(env.get("ty_top",0.0)))

    fixed=[]
    if "supports" in env:
        for n,fx,fy in np.array(env["supports"], dtype=int):
            if fx: fixed.append(2*int(n))
            if fy: fixed.append(2*int(n)+1)
    edges={"clamp_left":left,"clamp_right":right,"clamp_bottom":bottom,"clamp_top":top}
    for key,ids in edges.items():
        if int(env.get(key,0))==1:
            for i in ids: fixed.extend([2*i,2*i+1])
    fixed=np.array(sorted(set(fixed)),dtype=int)
    ndof=2*nn; free=np.setdiff1d(np.arange(ndof), fixed)
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
        eps=B@U[dofs]; sig=D@eps
        stress.append(sig)
    return {"nodes":nodes,"elems":elems,"U":U,"stress":np.array(stress)}
