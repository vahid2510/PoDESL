from __future__ import annotations
import numpy as np
from .vtk_writer import write_vtk_quad2d

def _von_mises_ps(sig):
    arr = np.asarray(sig, dtype=float)
    sx = arr[..., 0]
    sy = arr[..., 1]
    txy = arr[..., 2]
    return np.sqrt(sx * sx - sx * sy + sy * sy + 3.0 * txy * txy)

def _nodal_average_from_gp(elems, values, nn, n_gp_elem):
    accum = np.zeros((nn, values.shape[1]), dtype=float)
    counts = np.zeros(nn, dtype=float)
    for ee, conn in enumerate(elems):
        start = ee * n_gp_elem
        block = values[start:start + n_gp_elem]
        avg = block.mean(axis=0)
        for nid in conn:
            accum[nid] += avg
            counts[nid] += 1.0
    counts[counts == 0.0] = 1.0
    return accum / counts[:, None]
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
def _D_elastic(E, nu):
    c = E/(1-nu**2)
    return c*np.array([[1, nu, 0],[nu, 1, 0],[0, 0, (1-nu)/2]], dtype=float)
def _vm_eq_stress(sig):
    sx, sy, txy = sig
    return np.sqrt(sx*sx - sx*sy + sy*sy + 3*txy*txy)
def _G(E,nu): return E/(2*(1+nu))
def _ep_tangent(E, nu, H):
    G = _G(E,nu)
    Et = (3*G*H)/(3*G + H + 1e-30)
    c = E/(1-nu**2)
    De = c*np.array([[1, nu, 0],[nu, 1, 0],[0, 0, (1-nu)/2]], dtype=float)
    D = De.copy()
    D[2,2] = Et/(2*(1+nu) + 1e-30)
    r = max(0.0, min(1.0, Et/E))
    D[0,0] = De[0,0]* (0.5+0.5*r)
    D[1,1] = De[1,1]* (0.5+0.5*r)
    D[0,1] = De[0,1]* (0.5+0.5*r); D[1,0] = D[0,1]
    return D
def solve_ps_vonmises_fe(env):
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
    sigy0=float(env.get("sigma_y0", 250e6)); H=float(env.get("H", 1.0e9))
    ux_bar = env.get("ux_bar", None)
    traction_right_x = float(env.get("traction_right_x", 0.0))
    traction_right_y = float(env.get("traction_right_y", 0.0))
    steps=int(env.get("load_steps", 50)); max_iter=int(env.get("max_iter", 50))
    tol=float(env.get("tol", 1e-9)); penalty=float(env.get("penalty", 1e7))*E
    nn=nodes.shape[0]; ndof=2*nn
    left_nodes = [k for k,(x,y) in enumerate(nodes) if abs(x-0.0)<1e-12]
    right_nodes = [k for k,(x,y) in enumerate(nodes) if abs(x-Lx)<1e-12]
    fixed = []
    for n in left_nodes: fixed.extend([2*n, 2*n+1])
    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed)
    qp, qw = _gauss2(); n_gp_elem = len(qp)
    n_gp = len(elems)*n_gp_elem
    sig = np.zeros((n_gp,3)); alpha = np.zeros(n_gp); ep = np.zeros(n_gp)
    lam = 0.0; dlam = 1.0/steps; path = []
    def assemble(U, lam):
        K = np.zeros((ndof,ndof)); R = np.zeros(ndof)
        gp_idx = 0
        for e,conn in enumerate(elems):
            xe=nodes[conn,:]
            Ke=np.zeros((8,8)); Re=np.zeros(8)
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
                dofs=[]; 
                for n in conn: dofs.extend([2*n,2*n+1])
                ue = U[dofs]
                eps = B@ue
                De = _D_elastic(E,nu)
                sig_trial = De@eps
                seq = _vm_eq_stress(sig_trial)
                f = seq - (sigy0 + H*alpha[gp_idx])
                if f <= 0.0:
                    D = De; sig_gp = sig_trial
                else:
                    D = _ep_tangent(E,nu,H)
                    sig_gp = D@eps
                    dgamma = f/(3*_G(E,nu)+H + 1e-30)
                    alpha[gp_idx] += dgamma; ep[gp_idx] += dgamma
                sig[gp_idx,:] = sig_gp
                Ke += (B.T@D@B)*detJ*wt*t
                Re += (B.T@sig_gp)*detJ*wt*t
                gp_idx += 1
            for i in range(8):
                I=dofs[i]; R[I]+=Re[i]
                for j in range(8):
                    J=dofs[j]; K[I,J]+=Ke[i,j]
        if ux_bar is not None:
            ux_target = float(ux_bar)
            for n in right_nodes:
                I = 2*n
                R[I] += - penalty*(U[I] - lam*ux_target)
                K[I,I] += penalty
        if abs(traction_right_x) + abs(traction_right_y) > 0.0:
            idxs = sorted(right_nodes, key=lambda k: nodes[k,1])
            for a,b in zip(idxs[:-1], idxs[1:]):
                Le = np.linalg.norm(nodes[b]-nodes[a])
                fx = traction_right_x * t * Le / 2.0
                fy = traction_right_y * t * Le / 2.0
                R[2*a]   += -fx; R[2*a+1] += -fy
                R[2*b]   += -fx; R[2*b+1] += -fy
        return K, R
    U = np.zeros(ndof)
    for step in range(1, steps+1):
        target = min(1.0, lam + dlam)
        res = 0.0
        for it in range(max_iter):
            K, R = assemble(U, target)
            res = np.linalg.norm(R[free], ord=np.inf)
            if res < tol: break
            dU = np.zeros_like(U)
            dU[free] = np.linalg.solve(K[free][:,free], -R[free])
            U += dU
        lam = target
        path.append([step, lam, float(res)])
    # Post: Von Mises at GP and nodes
    vm_gp = _von_mises_ps(sig)
    node_sigma = _nodal_average_from_gp(elems, sig, nn, n_gp_elem)
    node_ep    = _nodal_average_from_gp(elems, ep.reshape(-1,1), nn, n_gp_elem).reshape(-1)
    node_vm    = _von_mises_ps(node_sigma)
    out = {"U": U, "path": np.array(path), "sigma_gp": sig, "alpha_gp": alpha, "ep_gp": ep,
           "sigma": node_sigma, "sigma_vm": node_vm, "ep": node_ep, "nodes": nodes, "elems": elems,
           "sigma_vm_gp": vm_gp}
    return out
