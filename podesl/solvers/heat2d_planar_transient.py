from __future__ import annotations
import numpy as np

def _q4_shape(xi, eta):
    N = 0.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)], float)
    dN_dxi = 0.25*np.array([[-(1-eta),-(1-xi)],[+(1-eta),-(1+xi)],[+(1+eta),+(1+xi)],[-(1+eta),+(1-xi)]], float)
    return N, dN_dxi
def solve_heat2d_planar_transient(env):
    nodes=np.asarray(env.get('nodes'), float); elems=np.asarray(env.get('elems'), int)
    k=float(env.get('k')); rho=float(env.get('rho')); cp=float(env.get('cp'))
    dt=float(env.get('dt')); t_end=float(env.get('t_end')); theta=float(env.get('theta',1.0))
    q=float(env.get('q',0.0)); nn=nodes.shape[0]
    K=np.zeros((nn,nn)); M=np.zeros((nn,nn)); F=np.zeros(nn)
    gp=[-1/np.sqrt(3), 1/np.sqrt(3)]
    for conn in elems:
        xe=nodes[conn,:]; Ke=np.zeros((4,4)); Me=np.zeros((4,4)); fe=np.zeros(4)
        for xi in gp:
            for eta in gp:
                N,dN_dxi=_q4_shape(xi,eta)
                J=np.zeros((2,2))
                for a in range(4):
                    J[0,0]+=dN_dxi[a,0]*xe[a,0]; J[0,1]+=dN_dxi[a,0]*xe[a,1]
                    J[1,0]+=dN_dxi[a,1]*xe[a,0]; J[1,1]+=dN_dxi[a,1]*xe[a,1]
                detJ=np.linalg.det(J); 
                if detJ<=0: raise ValueError('Invalid Q4 mapping')
                invJ=np.linalg.inv(J); dN_dx=dN_dxi @ invJ.T
                w=detJ; B=dN_dx
                Ke += k*(B @ B.T)*w
                Me += rho*cp*np.outer(N,N)*w
                fe += q*N*w
        for i,I in enumerate(conn):
            F[I]+=fe[i]
            for j,J in enumerate(conn):
                K[I,J]+=Ke[i,j]; M[I,J]+=Me[i,j]
    for spec in env.get('hedge', []):
        ee=int(spec[0]); edge=int(spec[1]); h=float(spec[2]); Tinf=float(spec[3])
        conn=elems[ee]; a,b=[(0,1),(1,2),(2,3),(3,0)][edge]
        n1,n2=conn[a],conn[b]; x1,y1=nodes[n1]; x2,y2=nodes[n2]
        L=float(np.hypot(x2-x1,y2-y1))
        Ke=(h*L/6.0)*np.array([[2,1],[1,2]],float); Fe=(h*Tinf*L/2.0)*np.array([1.0,1.0],float)
        idx=[n1,n2]
        for i,I in enumerate(idx):
            F[I]+=Fe[i]
            for j,J in enumerate(idx):
                K[I,J]+=Ke[i,j]
    fixed={}
    for ent in env.get('fix', []):
        n=int(ent[0]); val=float(ent[1]); fixed[n]=val
    fixed_idx=np.array(sorted(fixed.keys()),dtype=int); free_idx=np.setdiff1d(np.arange(nn),fixed_idx)
    T0=env.get('T0', None)
    Tn=(float(env.get('Tinit',300.0))*np.ones(nn)) if T0 is None else np.asarray(T0,float).reshape(nn)
    nsteps=int(np.ceil(t_end/dt)); times=np.linspace(0.0, nsteps*dt, nsteps+1); Thist=np.zeros((nsteps+1,nn)); Thist[0,:]=Tn
    A=M + theta*dt*K; B=M - (1-theta)*dt*K
    if fixed_idx.size>0:
        Af=A[np.ix_(free_idx,free_idx)]; Acf=A[np.ix_(free_idx,fixed_idx)]
    else:
        Af=A; Acf=None
    for s in range(1,nsteps+1):
        rhs=(B @ Tn) + dt*F
        if fixed_idx.size>0:
            rhs_free=rhs[free_idx] - Acf @ np.array([fixed[i] for i in fixed_idx])
            Tfree=np.linalg.solve(Af, rhs_free)
            Tn1=np.zeros(nn); Tn1[free_idx]=Tfree
            for i in fixed_idx: Tn1[i]=fixed[i]
        else:
            Tn1=np.linalg.solve(A, rhs)
        Thist[s,:]=Tn1; Tn=Tn1
    return {'T_hist':Thist,'times':times,'nodes':nodes}
