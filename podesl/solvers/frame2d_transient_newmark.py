import numpy as np
from ..linalg_utils import safe_solve
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
def _m_local(rho,A,I,L):
    m = rho*A*L/420.0 * np.array([
        [140,    0,          0,        70,     0,          0],
        [  0,  156,      22*L,         0,    54,     -13*L],
        [  0, 22*L,    4*L*L,          0, 13*L,    -3*L*L],
        [ 70,    0,          0,       140,     0,          0],
        [  0,   54,      13*L,         0,   156,     -22*L],
        [  0, -13*L,   -3*L*L,         0, -22*L,     4*L*L],
    ], dtype=float)
    return m
def solve_frame2d_transient_newmark(env):
    nodes = np.array(env["nodes"], dtype=float)
    beams = np.array(env["elems"], dtype=int)
    E=float(env["E"]); A=float(env["A"]); I=float(env["I"]); rho=float(env["rho"])
    dt=float(env.get("dt", 0.002)); t_end=float(env.get("t_end", 1.0))
    beta=float(env.get("beta", 0.25)); gamma=float(env.get("gamma", 0.5))
    a0=float(env.get("rayleigh_a0", 0.0)); a1=float(env.get("rayleigh_a1", 0.0))
    nn=nodes.shape[0]; ndof=3*nn
    K=np.zeros((ndof,ndof)); M=np.zeros((ndof,ndof))
    for (i,j) in beams:
        xi, yi = nodes[i]; xj, yj = nodes[j]
        dx=xj-xi; dy=yj-yi; L=float(np.hypot(dx,dy)); c=dx/L; s=dy/L
        T=_T(c,s)
        k_loc=_k_local(E,A,I,L); m_loc=_m_local(rho,A,I,L)
        k_gl=T.T@k_loc@T; m_gl=T.T@m_loc@T
        dofs=[3*i,3*i+1,3*i+2, 3*j,3*j+1,3*j+2]
        for p in range(6):
            for q in range(6):
                K[dofs[p],dofs[q]]+=k_gl[p,q]
                M[dofs[p],dofs[q]]+=m_gl[p,q]
    C = a0*M + a1*K
    nsteps=int(np.ceil(t_end/dt))
    times=np.linspace(0,nsteps*dt,nsteps+1)
    entries = np.array(env.get("loads_t", []), dtype=float).reshape(-1,5) if len(env.get("loads_t", []))>0 else np.zeros((0,5))
    hist = {}
    for row in entries:
        t, n, fx, fy, mz = row
        n=int(n)
        for dof,val in [(3*n,fx),(3*n+1,fy),(3*n+2,mz)]:
            if dof not in hist: hist[dof]=[]
            hist[dof].append((t,val))
    for dof in list(hist.keys()):
        hist[dof]=sorted(hist[dof], key=lambda x:x[0])
    def F_at(t):
        F=np.zeros(ndof)
        for dof, tv in hist.items():
            if t<=tv[0][0]: F[dof]=tv[0][1]; continue
            if t>=tv[-1][0]: F[dof]=tv[-1][1]; continue
            for k in range(len(tv)-1):
                t0,v0=tv[k]; t1,v1=tv[k+1]
                if t0<=t<=t1:
                    w=(t-t0)/(t1-t0); F[dof]=v0*(1-w)+v1*w; break
        return F
    fixed=[]
    for n,fu,fv,fr in env.get("supports", []):
        n=int(n)
        if fu: fixed.append(3*n)
        if fv: fixed.append(3*n+1)
        if fr: fixed.append(3*n+2)
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check supports.")
    U=np.zeros(ndof); V=np.zeros(ndof); A=np.zeros(ndof)
    Keff = K + gamma/(beta*dt)*C + M/(beta*dt*dt)
    U_hist=[]; V_hist=[]; A_hist=[]
    for it,t in enumerate(times):
        F = F_at(t)
        Feff = F              + M@( (1/(beta*dt*dt))*U + (1/(beta*dt))*V + (1/(2*beta)-1)*A )              + C@( (gamma/(beta*dt))*U + (gamma/beta -1)*V + dt*(gamma/(2*beta)-1)*A )
        U_new = U.copy()
        U_new[free] = safe_solve(Keff[free][:,free], Feff[free])
        A_new = (1/(beta*dt*dt))*(U_new - U) - (1/(beta*dt))*V - ((1/(2*beta))-1)*A
        V_new = V + dt*((1-gamma)*A + gamma*A_new)
        U,V,A = U_new, V_new, A_new
        U_hist.append(U.copy()); V_hist.append(V.copy()); A_hist.append(A.copy())
    return {"times":times, "U_hist":np.array(U_hist), "V_hist":np.array(V_hist), "A_hist":np.array(A_hist)}
