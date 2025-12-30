import numpy as np
from .nl_core import NLControls, solve_load_steps, solve_arc_length

def _T(c,s):
    R = np.array([[c, s, 0],
                  [-s, c, 0],
                  [0, 0, 1]], dtype=float)
    T = np.zeros((6,6)); T[:3,:3]=R; T[3:,3:]=R
    return T
def _k_local(E,A,I,L):
    k = np.array([
        [ A*E/L,       0,           0,      -A*E/L,       0,           0],
        [ 0,     12*E*I/L**3,  6*E*I/L**2,   0,   -12*E*I/L**3,  6*E*I/L**2],
        [ 0,      6*E*I/L**2,  4*E*I/L,      0,   -6*E*I/L**2,  2*E*I/L   ],
        [-A*E/L,       0,           0,       A*E/L,       0,           0],
        [ 0,    -12*E*I/L**3, -6*E*I/L**2,  0,    12*E*I/L**3, -6*E*I/L**2],
        [ 0,      6*E*I/L**2,  2*E*I/L,     0,   -6*E*I/L**2,  4*E*I/L   ]
    ], dtype=float)
    return k
def _kg_local(N,L):
    c1 = N/(30*L)
    Kg = c1 * np.array([
        [ 0,   0,    0,   0,   0,    0],
        [ 0,  36,   3*L,  0, -36,   3*L],
        [ 0,  3*L,  4*L*L,0, -3*L, -L*L],
        [ 0,   0,    0,   0,   0,    0],
        [ 0, -36,  -3*L,  0,  36,  -3*L],
        [ 0,  3*L, -L*L,  0, -3*L,  4*L*L]
    ], dtype=float)
    return Kg

def solve_frame2d_static_nl(env):
    nodes = np.array(env["nodes"], dtype=float)
    beams = np.array(env["elems"], dtype=int)
    E = float(env["E"]); A = float(env["A"]); I = float(env["I"])
    nn = nodes.shape[0]; ndof=3*nn
    Fext_full = np.zeros(ndof)
    for n,fx,fy,mz in env.get("loads", []):
        n=int(n); Fext_full[3*n]+=float(fx); Fext_full[3*n+1]+=float(fy); Fext_full[3*n+2]+=float(mz)
    fixed=[]
    for n,fu,fv,fr in env.get("supports", []):
        n=int(n)
        if fu: fixed.append(3*n)
        if fv: fixed.append(3*n+1)
        if fr: fixed.append(3*n+2)
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")

    Cmap = env.get("CONTROLS", {})
    controls = NLControls(
        method=str(Cmap.get("method", "loadstep")).lower(),
        load_steps=int(Cmap.get("load_steps", env.get("load_steps", 20))),
        max_iter=int(Cmap.get("max_iter", env.get("max_iter", 40))),
        tol=float(Cmap.get("tol", env.get("tol", 1e-9))),
        line_search=bool(Cmap.get("line_search", True)),
        arc_len=float(Cmap.get("arc_len", 1e-3)),
        arc_minmax=tuple(Cmap.get("arc_minmax", [1e-5, 1e-2])),
        report_each=int(Cmap.get("report_each", 1)),
        sparse=bool(env.get("SPARSE", Cmap.get("sparse", False)))
    )

    def assemble(U):
        ndof = U.size
        K=np.zeros((ndof,ndof)); Fint=np.zeros(ndof)
        for (i,j) in beams:
            xi, yi = nodes[i]; xj, yj = nodes[j]
            dx=xj-xi; dy=yj-yi; L=float(np.hypot(dx,dy)); c=dx/L; s=dy/L
            T=_T(c,s)
            dofs=[3*i,3*i+1,3*i+2, 3*j,3*j+1,3*j+2]
            ue = U[dofs]; ul = T@ue
            N = E*A*(ul[3]-ul[0])/L
            k_loc=_k_local(E,A,I,L); kg_loc=_kg_local(N,L)
            k_gl=T.T@k_loc@T + T.T@kg_loc@T
            fint = (T.T@k_loc@T)@ue
            for p in range(6): Fint[dofs[p]]+=fint[p]
            for p in range(6):
                for q in range(6):
                    K[dofs[p],dofs[q]]+=k_gl[p,q]
        return K, Fint

    if controls.method == "arc":
        U, path, _ = solve_arc_length(assemble, ndof, free, Fext_full, controls)
    else:
        U, path, _ = solve_load_steps(assemble, ndof, free, Fext_full, controls)
    K, Fint = assemble(U)
    Reac = (K@U - Fext_full)
    return {"U":U, "reactions":Reac, "path":path}
