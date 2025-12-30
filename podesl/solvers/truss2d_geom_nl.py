import numpy as np
from .nl_core import NLControls, solve_load_steps, solve_arc_length

def solve_truss2d_geom_nl(env):
    nodes = np.array(env["nodes"], dtype=float)
    bars  = np.array(env["elems"], dtype=int)
    E = float(env["E"]); A = float(env["A"])
    nn = nodes.shape[0]; ndof=2*nn
    Fext_full = np.zeros(ndof)
    for n,fx,fy in env.get("loads", []):
        n=int(n); Fext_full[2*n]+=float(fx); Fext_full[2*n+1]+=float(fy)
    fixed=[]
    for n,fx,fy in env.get("supports", []):
        n=int(n)
        if fx: fixed.append(2*n)
        if fy: fixed.append(2*n+1)
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    L0 = np.array([np.hypot(*(nodes[j]-nodes[i])) for (i,j) in bars], dtype=float)

    Cmap = env.get("CONTROLS", {})
    controls = NLControls(
        method=str(Cmap.get("method", "loadstep")).lower(),
        load_steps=int(Cmap.get("load_steps", env.get("load_steps", 20))),
        max_iter=int(Cmap.get("max_iter", env.get("max_iter", 40))),
        tol=float(Cmap.get("tol", env.get("tol", 1e-10))),
        line_search=bool(Cmap.get("line_search", True)),
        arc_len=float(Cmap.get("arc_len", 1e-3)),
        arc_minmax=tuple(Cmap.get("arc_minmax", [1e-5, 1e-2])),
        report_each=int(Cmap.get("report_each", 1)),
        sparse=bool(env.get("SPARSE", Cmap.get("sparse", False)))
    )

    def assemble(U):
        ndof = U.size
        K=np.zeros((ndof,ndof)); Fint=np.zeros(ndof)
        for e,(i,j) in enumerate(bars):
            xi=nodes[i]; xj=nodes[j]
            ui=U[2*i:2*i+2]; uj=U[2*j:2*j+2]
            xi2=xi+ui; xj2=xj+uj
            dx=xj2-xi2; L=float(np.hypot(dx[0],dx[1])); t=dx/L
            eps=(L-L0[e])/L0[e]
            N=E*A*eps
            B=(1.0/L0[e])*np.array([-t[0],-t[1], t[0], t[1]])
            k_mat=(E*A)*np.outer(B,B)
            c,s=t[0],t[1]
            G = np.array([[ s*s, -s*c, -s*s,  s*c],
                          [-s*c,  c*c,  s*c, -c*c],
                          [-s*s,  s*c,  s*s, -s*c],
                          [ s*c, -c*c, -s*c,  c*c]])
            k_geo=(N/L)*G
            fint=(N/L)*np.array([-dx[0],-dx[1], dx[0], dx[1]])
            dofs=[2*i,2*i+1,2*j,2*j+1]
            for p in range(4): Fint[dofs[p]]+=fint[p]
            ke=k_mat+k_geo
            for p in range(4):
                for q in range(4):
                    K[dofs[p],dofs[q]]+=ke[p,q]
        return K, Fint

    if controls.method == "arc":
        U, path, _ = solve_arc_length(assemble, ndof, free, Fext_full, controls)
    else:
        U, path, _ = solve_load_steps(assemble, ndof, free, Fext_full, controls)

    K, Fint = assemble(U)
    Reac = (K@U - Fext_full)
    return {"U":U, "reactions":Reac, "path":path}
