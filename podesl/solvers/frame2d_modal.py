import numpy as np
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
def _m_local(rho,A,I,L):
    m = rho*A*L
    M = np.array([
        [ m/3,     0,         0,     m/6,     0,         0],
        [ 0,  13*m/35,  11*m*L/210,   0,   9*m/70, -13*m*L/420],
        [ 0, 11*m*L/210,  m*L**2/105, 0,  13*m*L/420, -m*L**2/140],
        [ m/6,    0,         0,     m/3,      0,         0],
        [ 0,   9*m/70,  13*m*L/420,  0,  13*m/35, -11*m*L/210],
        [ 0,-13*m*L/420, -m*L**2/140,0,-11*m*L/210,  m*L**2/105]
    ], dtype=float)
    return M
def _T(c,s):
    R = np.array([[c, s, 0],
                  [-s, c, 0],
                  [0, 0, 1]], dtype=float)
    T = np.zeros((6,6)); T[:3,:3]=R; T[3:,3:]=R
    return T
def solve_frame2d_modal(env):
    nodes = np.array(env["nodes"], dtype=float)
    beams = np.array(env["elems"], dtype=int)
    E = float(env["E"]); A = float(env["A"]); I = float(env["I"]); rho=float(env["rho"])
    nn = nodes.shape[0]; ndof=3*nn
    K=np.zeros((ndof,ndof)); M=np.zeros((ndof,ndof))
    for e,(i,j) in enumerate(beams):
        xi, yi = nodes[i]; xj, yj = nodes[j]
        dx=xj-xi; dy=yj-yi; L=float(np.hypot(dx,dy)); c=dx/L; s=dy/L
        k_loc=_k_local(E,A,I,L); m_loc=_m_local(rho,A,I,L); T=_T(c,s)
        k_gl=T.T@k_loc@T; m_gl=T.T@m_loc@T
        dofs=[3*i,3*i+1,3*i+2, 3*j,3*j+1,3*j+2]
        for p in range(6):
            for q in range(6):
                K[dofs[p],dofs[q]]+=k_gl[p,q]
                M[dofs[p],dofs[q]]+=m_gl[p,q]
    fixed=[]
    for n,fu,fv,fr in env.get("supports", []):
        n=int(n)
        if fu: fixed.append(3*n)
        if fv: fixed.append(3*n+1)
        if fr: fixed.append(3*n+2)
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    MinvK = np.linalg.solve(M[free][:,free], K[free][:,free])
    vals, vecs = np.linalg.eig(MinvK)
    vals=np.real(vals); vecs=np.real(vecs)
    idx=np.argsort(vals)
    w2=vals[idx]; Phi=vecs[:,idx]
    w=np.sqrt(np.clip(w2,0,None)); f=w/(2*np.pi)
    Phi_full=np.zeros((ndof, Phi.shape[1])); Phi_full[free,:]=Phi
    return {"omega":w, "freq":f, "modes":Phi_full}
