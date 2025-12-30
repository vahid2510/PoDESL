import numpy as np
def _elem_k(E,A,x1,y1,x2,y2):
    dx=x2-x1; dy=y2-y1; L=np.hypot(dx,dy); c=dx/L; s=dy/L
    k = (E*A/L)*np.array([[ c*c,  c*s, -c*c, -c*s],
                          [ c*s,  s*s, -c*s, -s*s],
                          [-c*c, -c*s,  c*c,  c*s],
                          [-c*s, -s*s,  c*s,  s*s]], dtype=float)
    return k
def solve_truss2d_static(env):
    nodes = np.array(env["nodes"], dtype=float)  # [[x,y],...]
    bars  = np.array(env["elems"], dtype=int)    # [[n1,n2],...]
    E = float(env["E"])
    A = env.get("A", None)
    if isinstance(A, (list, tuple)) and len(A)==len(bars):
        A = np.array(A, dtype=float)
    else:
        A = float(A if A is not None else env.get("A_default", 1.0))
    nn = nodes.shape[0]; ndof=2*nn
    K = np.zeros((ndof,ndof)); F = np.zeros(ndof)
    for e,(i,j) in enumerate(bars):
        a = A[e] if isinstance(A, np.ndarray) else A
        ke = _elem_k(E,a,*nodes[i],*nodes[j])
        dofs=[2*i,2*i+1,2*j,2*j+1]
        for p in range(4):
            for q in range(4):
                K[dofs[p],dofs[q]] += ke[p,q]
    if "loads" in env:
        loads = np.array(env["loads"], dtype=float)
        for n,Fx,Fy in loads:
            n=int(n); F[2*n]+=Fx; F[2*n+1]+=Fy
    fixed=[]
    if "supports" in env:
        sp = np.array(env["supports"], dtype=int)
        for n,fx,fy in sp:
            if fx: fixed.append(2*int(n))
            if fy: fixed.append(2*int(n)+1)
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof),fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    U=np.zeros(ndof)
    U[free]=np.linalg.solve(K[free][:,free],F[free])
    # axial forces
    N=[]
    for e,(i,j) in enumerate(bars):
        ui=U[2*i:2*i+2]; uj=U[2*j:2*j+2]
        xi=nodes[i]; xj=nodes[j]
        dx=xj[0]-xi[0]; dy=xj[1]-xi[1]; L=np.hypot(dx,dy); c=dx/L; s=dy/L
        du = (uj-ui) @ np.array([c,s])
        a = A[e] if isinstance(A, np.ndarray) else A
        N.append(E*a*du/L)
    return {"U":U, "N":np.array(N)}
