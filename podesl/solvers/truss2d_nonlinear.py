import numpy as np
from ..linalg_utils import safe_solve
def _geom_terms(xi,xj,ui,uj,E,A,L0):
    xi2 = xi + ui; xj2 = xj + uj
    dx = xj2 - xi2; L = float(np.hypot(dx[0], dx[1])); t = dx / L
    eps = (L - L0)/L0
    N = E*A*eps
    B = (1.0/L0) * np.array([-t[0],-t[1], t[0], t[1]])
    k_mat = E*A * np.outer(B, B)
    c,s = t[0], t[1]
    c2, s2, cs = c*c, s*s, c*s
    G = np.array([[ s2, -cs, -s2,  cs],
                  [-cs,  c2,  cs, -c2],
                  [-s2,  cs,  s2, -cs],
                  [ cs, -c2, -cs,  c2]])
    k_geo = (N / L) * G
    f_int = (N/L) * np.array([-dx[0], -dx[1], dx[0], dx[1]])
    return f_int, (k_mat + k_geo)

def solve_truss2d_nonlinear(env):
    nodes = np.array(env["nodes"], dtype=float)
    bars  = np.array(env["elems"], dtype=int)
    E = float(env["E"])
    A = env.get("A", 1.0)
    if isinstance(A, (list, tuple)) and len(A)==len(bars):
        A = np.array(A, dtype=float)
    else:
        A = float(A)
    nn = nodes.shape[0]; ndof=2*nn
    U=np.zeros(ndof)
    Fext = np.zeros(ndof)
    if "loads" in env:
        for n,Fx,Fy in np.array(env["loads"], dtype=float):
            n=int(n); Fext[2*n]+=Fx; Fext[2*n+1]+=Fy
    fixed=[]
    if "supports" in env:
        for n,fx,fy in np.array(env["supports"], dtype=int):
            if fx: fixed.append(2*int(n))
            if fy: fixed.append(2*int(n)+1)
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof),fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    L0 = np.array([np.hypot(*(nodes[j]-nodes[i])) for (i,j) in bars], dtype=float)

    max_iter = int(env.get("max_iter", 50))
    tol = float(env.get("tol", 1e-8))
    for it in range(max_iter):
        K = np.zeros((ndof,ndof)); Fint = np.zeros(ndof)
        for e,(i,j) in enumerate(bars):
            ui = U[2*i:2*i+2]; uj = U[2*j:2*j+2]
            xi = nodes[i]; xj = nodes[j]
            a = A[e] if isinstance(A, np.ndarray) else A
            f_e, k_e = _geom_terms(xi,xj,ui,uj,E,a,L0[e])
            dofs=[2*i,2*i+1,2*j,2*j+1]
            for p in range(4): Fint[dofs[p]] += f_e[p]
            for p in range(4):
                for q in range(4):
                    K[dofs[p],dofs[q]] += k_e[p,q]
        R = Fext - Fint
        res = np.linalg.norm(R[free], ord=np.inf)
        if res < tol: break
        du = np.zeros(ndof)
        du[free] = safe_solve(K[free][:,free], R[free])
        U += du
        if np.linalg.norm(du[free], ord=np.inf) < tol: break

    N=[]
    for e,(i,j) in enumerate(bars):
        xi = nodes[i]; xj = nodes[j]
        ui = U[2*i:2*i+2]; uj = U[2*j:2*j+2]
        L = np.hypot(*(xi+ui - (xj+uj)))
        eps = (L - L0[e])/L0[e]
        a = A[e] if isinstance(A, np.ndarray) else A
        N.append(E*a*eps)
    return {"U":U, "N":np.array(N), "iters": it+1}
