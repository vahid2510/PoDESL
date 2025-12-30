
from __future__ import annotations
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

def _D_plane_strain(E, nu):
    lam = E*nu/((1+nu)*(1-2*nu))
    mu  = E/(2*(1+nu))
    return np.array([[lam+2*mu, lam,       0.0],
                     [lam,       lam+2*mu, 0.0],
                     [0.0,       0.0,      mu ]], dtype=float)

def _von_mises_ps(sig):
    sx, sy, txy = sig
    return np.sqrt(sx*sx - sx*sy + sy*sy + 3*txy*txy)

def solve_ps_linear_fe(env):
    Lx=float(env.get("Lx",1.0)); Ly=float(env.get("Ly",1.0)); t=float(env.get("t",1.0))
    nx=int(env.get("nx",8)); ny=int(env.get("ny",8))
    xs = np.linspace(0.0, Lx, nx+1); ys = np.linspace(0.0, Ly, ny+1)
    nodes = np.array([[x,y] for j in range(ny+1) for x in xs for y in [ys[j]]], dtype=float)
    def nid(i,j): return j*(nx+1) + i
    elems=[]
    for j in range(ny):
        for i in range(nx):
            n00=nid(i,j); n10=nid(i+1,j); n11=nid(i+1,j+1); n01=nid(i,j+1)
            elems.append([n00,n10,n11,n01])
    elems=np.array(elems,dtype=int)

    E=float(env.get("E",210e9)); nu=float(env.get("nu",0.3))
    D = _D_plane_strain(E, nu)

    traction_right_x = float(env.get("traction_right_x", 0.0))
    traction_right_y = float(env.get("traction_right_y", 0.0))
    clamp_left = bool(env.get("clamp_left", True))

    nn=nodes.shape[0]; ndof=2*nn
    K = np.zeros((ndof,ndof)); F = np.zeros(ndof)

    qp, qw = _gauss2()

    for e,conn in enumerate(elems):
        xe=nodes[conn,:]
        Ke=np.zeros((8,8))
        for (xi,eta),wt in zip(qp,qw):
            N, dN_dxi = _shape_Q4(xi,eta)
            J = dN_dxi@xe; detJ=np.linalg.det(J); invJ=np.linalg.inv(J)
            dN_dx = invJ@dN_dxi
            B=np.zeros((3,8))
            for i in range(4):
                B[0,2*i]   = dN_dx[0,i]
                B[1,2*i+1] = dN_dx[1,i]
                B[2,2*i]   = dN_dx[1,i]
                B[2,2*i+1] = dN_dx[0,i]
            Ke += (B.T@D@B)*detJ*wt*t
        dofs=[]; 
        for n in conn: dofs.extend([2*n,2*n+1])
        for i in range(8):
            I=dofs[i]
            for j in range(8):
                J=dofs[j]; K[I,J]+=Ke[i,j]

    right_nodes = [k for k,(x,y) in enumerate(nodes) if abs(x-Lx)<1e-12]
    if abs(traction_right_x)+abs(traction_right_y)>0.0:
        idxs = sorted(right_nodes, key=lambda k: nodes[k,1])
        for a,b in zip(idxs[:-1], idxs[1:]):
            Le = np.linalg.norm(nodes[b]-nodes[a])
            fx = traction_right_x * t * Le / 2.0
            fy = traction_right_y * t * Le / 2.0
            F[2*a]   += fx;  F[2*a+1] += fy
            F[2*b]   += fx;  F[2*b+1] += fy

    fixed = []
    if clamp_left:
        for k,(x,y) in enumerate(nodes):
            if abs(x-0.0)<1e-12:
                fixed.extend([2*k,2*k+1])
    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed)

    U = np.zeros(ndof)
    if free.size>0:
        U[free] = np.linalg.solve(K[free][:,free], F[free])

    sig_gp = np.zeros((len(elems),3)); vm_gp = np.zeros(len(elems))
    for e,conn in enumerate(elems):
        xe = nodes[conn,:]
        xi=eta=0.0
        N, dN_dxi = _shape_Q4(xi,eta)
        J = dN_dxi@xe; invJ=np.linalg.inv(J)
        dN_dx = invJ@dN_dxi
        B=np.zeros((3,8))
        for i in range(4):
            B[0,2*i]   = dN_dx[0,i]
            B[1,2*i+1] = dN_dx[1,i]
            B[2,2*i]   = dN_dx[1,i]
            B[2,2*i+1] = dN_dx[0,i]
        dofs=[]; 
        for n in conn: dofs.extend([2*n,2*n+1])
        ue = U[dofs]
        eps = B@ue
        sig = D@eps
        sig_gp[e,:] = sig
        vm_gp[e] = _von_mises_ps(sig)

    nn=nodes.shape[0]
    node_sigma = np.zeros((nn,3)); node_vm = np.zeros(nn); cnt = np.zeros(nn)
    for e,conn in enumerate(elems):
        for n in conn:
            node_sigma[n,:] += sig_gp[e,:]
            node_vm[n]      += vm_gp[e]
            cnt[n] += 1.0
    cnt[cnt==0]=1.0
    node_sigma /= cnt[:,None]
    node_vm    /= cnt

    return {"U": U, "F": F, "K": K, "sigma": node_sigma, "sigma_vm": node_vm,
            "nodes": nodes, "elems": elems}
