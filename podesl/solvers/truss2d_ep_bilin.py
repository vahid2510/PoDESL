import numpy as np
from ..linalg_utils import safe_solve
def solve_truss2d_ep_bilinear(env):
    nodes = np.array(env["nodes"], dtype=float)
    bars  = np.array(env["elems"], dtype=int)
    E = float(env["E"]); A = float(env["A"])
    sigy = float(env["sigy"]); H = float(env.get("H", 0.0))
    max_iter = int(env.get("max_iter", 60)); tol = float(env.get("tol", 1e-10))
    nn = nodes.shape[0]; ndof=2*nn
    U = np.zeros(ndof); Fext = np.zeros(ndof)
    for n,Fx,Fy in env.get("loads", []):
        n=int(n); Fext[2*n]+=float(Fx); Fext[2*n+1]+=float(Fy)
    fixed=[]
    for n,fx,fy in env.get("supports", []):
        if int(fx): fixed.append(2*int(n))
        if int(fy): fixed.append(2*int(n)+1)
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof),fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    ne = bars.shape[0]
    L0 = np.array([np.hypot(*(nodes[j]-nodes[i])) for (i,j) in bars], dtype=float)
    ep = np.zeros(ne); a_back = np.zeros(ne)
    def elem_tangent(i,j,ui,uj,e):
        xi=nodes[i]; xj=nodes[j]
        xi2=xi+ui; xj2=xj+uj
        dx=xj2 - xi2; L=float(np.hypot(dx[0],dx[1])); t=dx/L
        eps = (L - L0[e])/L0[e]
        sig_trial = E*(eps - ep[e])
        xi_v = sig_trial - a_back[e]
        ftrial = abs(xi_v) - sigy
        if ftrial<=0:
            N = sig_trial*A; Et = E
        else:
            sign = np.sign(xi_v) if xi_v!=0 else 1.0
            dgamma = ftrial/(E+H)
            ep[e] += dgamma*sign
            a_back[e] += H*dgamma*sign
            sig_trial = E*(eps - ep[e])
            N = (a_back[e] + sign*sigy)*A if H>0 else (sig_trial)*A
            Et = (E*H)/(E+H) if H>0 else 0.0
        B = (1.0/L0[e]) * np.array([-t[0],-t[1], t[0], t[1]])
        k_mat = (Et*A) * np.outer(B,B)
        c,s=t[0],t[1]
        G = np.array([[ s*s, -s*c, -s*s,  s*c],
                      [-s*c,  c*c,  s*c, -c*c],
                      [-s*s,  s*c,  s*s, -s*c],
                      [ s*c, -c*c, -s*c,  c*c]])
        k_geo = (N/L) * G
        fint = (N/L) * np.array([-dx[0], -dx[1], dx[0], dx[1]])
        return fint, k_mat + k_geo
    for it in range(max_iter):
        K = np.zeros((ndof,ndof)); Fint = np.zeros(ndof)
        for e,(i,j) in enumerate(bars):
            ui=U[2*i:2*i+2]; uj=U[2*j:2*j+2]
            fint, ke = elem_tangent(i,j,ui,uj,e)
            dofs=[2*i,2*i+1,2*j,2*j+1]
            for p in range(4): Fint[dofs[p]]+=fint[p]
            for p in range(4):
                for q in range(4):
                    K[dofs[p],dofs[q]]+=ke[p,q]
        R = Fext - Fint
        res = np.linalg.norm(R[free], ord=np.inf)
        if res<tol: break
        du = np.zeros(ndof); du[free]=safe_solve(K[free][:,free],R[free])
        U += du
        if np.linalg.norm(du[free], ord=np.inf)<tol: break
    # recompute axial forces
    N_list=[]
    for e,(i,j) in enumerate(bars):
        xi=nodes[i]; xj=nodes[j]
        ui=U[2*i:2*i+2]; uj=U[2*j:2*j+2]
        L = np.hypot(*(xj+uj - (xi+ui)))
        eps = (L - L0[e])/L0[e]
        sig = E*(eps - ep[e])
        N_list.append(sig*A)
    return {"U":U, "N":np.array(N_list), "ep":ep, "aback":a_back}
