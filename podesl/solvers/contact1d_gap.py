import numpy as np
def solve_contact1d_gap(env):
    nn = int(env.get("nn", 0))
    if "nodes" in env:
        nn = max(nn, len(env["nodes"]))
    if "springs" in env:
        for i,j,_ in env["springs"]:
            nn = max(nn, int(i)+1, int(j)+1)
    if "contact_gaps" in env:
        for i,j,_,_ in env["contact_gaps"]:
            nn = max(nn, int(i)+1, int(j)+1)
    ndof = nn
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    for i,j,k in env.get("springs", []):
        i=int(i); j=int(j); k=float(k)
        K[i,i]+=k; K[j,j]+=k; K[i,j]-=k; K[j,i]-=k

    for n,f in env.get("loads", []):
        F[int(n)] += float(f)

    fixed=[]
    for n,fx in env.get("supports", []):
        if int(fx): fixed.append(int(n))
    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed)

    U = np.zeros(ndof)
    if free.size>0:
        U[free] = np.linalg.solve(K[free][:,free], F[free])

    max_iter = int(env.get("max_iter", 50)); tol = float(env.get("tol", 1e-10))
    active = set()
    for _ in range(max_iter):
        changed = False
        for i,j,g0,kp in env.get("contact_gaps", []):
            i=int(i); j=int(j); g0=float(g0)
            key=(i,j)
            gap = (U[i]-U[j]) - g0
            if gap < 0 and key not in active:
                active.add(key); changed=True
            if gap >= 0 and key in active:
                active.remove(key); changed=True
        if not changed: break
        Kc = K.copy(); Fc = F.copy()
        for (i,j) in active:
            kp=None; g0=None
            for ii,jj,gg,kk in env.get("contact_gaps", []):
                if int(ii)==i and int(jj)==j:
                    g0=float(gg); kp=float(kk); break
            Kc[i,i]+=kp; Kc[j,j]+=kp; Kc[i,j]-=kp; Kc[j,i]-=kp
            Fc[i] += kp*g0; Fc[j] -= kp*g0
        U = np.zeros(ndof)
        if free.size>0:
            U[free] = np.linalg.solve(Kc[free][:,free], Fc[free])
        if np.linalg.norm((Kc@U - Fc)[free], ord=np.inf) < tol:
            break
    return {"U":U, "active_contacts": list(active)}
