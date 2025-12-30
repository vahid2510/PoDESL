from __future__ import annotations
import numpy as np

def solve_heat1d_transient(env):
    nodes=np.asarray(env.get('nodes'), float).reshape(-1)
    elems=np.asarray(env.get('elems'), int)
    k=float(env.get('k')); rho=float(env.get('rho')); cp=float(env.get('cp'))
    dt=float(env.get('dt')); t_end=float(env.get('t_end')); theta=float(env.get('theta',1.0))
    qvol=float(env.get('q', 0.0))
    nn=nodes.size
    K=np.zeros((nn,nn)); M=np.zeros((nn,nn)); F=np.zeros(nn)
    for (i,j) in elems:
        x1=nodes[i]; x2=nodes[j]; L=float(abs(x2-x1))
        Ke = k/L * np.array([[ 1,-1],[-1, 1]], float)
        Me = rho*cp*L/6.0 * np.array([[2,1],[1,2]], float)
        Fe = qvol*L/2.0 * np.array([1.0,1.0], float)
        idx=[i,j]
        for a,A in enumerate(idx):
            F[A]+=Fe[a]
            for b,B in enumerate(idx):
                K[A,B]+=Ke[a,b]; M[A,B]+=Me[a,b]
    for ent in env.get('hedge', []):
        ee=int(ent[0]); end=int(ent[1]); h=float(ent[2]); Tinf=float(ent[3])
        i,j = elems[ee]; L=float(abs(nodes[j]-nodes[i]))
        n = i if end==0 else j
        K[n,n]+=h; F[n]+=h*Tinf
    fixed={}
    for ent in env.get('fix', []):
        n=int(ent[0]); val=float(ent[1]); fixed[n]=val
    fixed_idx=np.array(sorted(fixed.keys()), dtype=int); free_idx=np.setdiff1d(np.arange(nn), fixed_idx)
    T0=env.get('T0', None)
    if T0 is None:
        Tn = float(env.get('Tinit',300.0))*np.ones(nn)
    else:
        T_arr = np.asarray(T0, float)
        if T_arr.size == 1:
            Tn = np.full(nn, float(T_arr))
        else:
            Tn = T_arr.reshape(nn)
    nsteps=int(np.ceil(t_end/dt)); times=np.linspace(0.0,nsteps*dt,nsteps+1); Thist=np.zeros((nsteps+1, nn)); Thist[0,:]=Tn
    A = M + theta*dt*K; B = M - (1-theta)*dt*K
    if fixed_idx.size>0:
        Af=A[np.ix_(free_idx, free_idx)]; Acf=A[np.ix_(free_idx, fixed_idx)]
    else:
        Af=A; Acf=None
    for s in range(1, nsteps+1):
        rhs=(B @ Tn) + dt*F
        if fixed_idx.size>0:
            rhs_f = rhs[free_idx] - Acf @ np.array([fixed[i] for i in fixed_idx])
            Tf = np.linalg.solve(Af, rhs_f)
            Tn1=np.zeros(nn); Tn1[free_idx]=Tf
            for i in fixed_idx: Tn1[i]=fixed[i]
        else:
            Tn1=np.linalg.solve(A, rhs)
        Thist[s,:]=Tn1; Tn=Tn1
    return {'T_hist':Thist, 'T': Thist, 'times':times, 'nodes':nodes, 'elems': env.get('elems')}
