from __future__ import annotations
import numpy as np

def _H8_shape(xi,eta,zeta):
    N = 0.125*np.array([(1-xi)*(1-eta)*(1-zeta),
                        (1+xi)*(1-eta)*(1-zeta),
                        (1+xi)*(1+eta)*(1-zeta),
                        (1-xi)*(1+eta)*(1-zeta),
                        (1-xi)*(1-eta)*(1+zeta),
                        (1+xi)*(1-eta)*(1+zeta),
                        (1+xi)*(1+eta)*(1+zeta),
                        (1-xi)*(1+eta)*(1+zeta)], float)
    dN = 0.125*np.array([[-(1-eta)*(1-zeta), -(1-xi)*(1-zeta), -(1-xi)*(1-eta)],
                         [ +(1-eta)*(1-zeta), -(1+xi)*(1-zeta), -(1+xi)*(1-eta)],
                         [ +(1+eta)*(1-zeta), +(1+xi)*(1-zeta), -(1+xi)*(1+eta)],
                         [ -(1+eta)*(1-zeta), +(1-xi)*(1-zeta), -(1-xi)*(1+eta)],
                         [ -(1-eta)*(1+zeta), -(1-xi)*(1+zeta), +(1-xi)*(1-eta)],
                         [ +(1-eta)*(1+zeta), -(1+xi)*(1+zeta), +(1+xi)*(1-eta)],
                         [ +(1+eta)*(1+zeta), +(1+xi)*(1+zeta), +(1+xi)*(1+eta)],
                         [ -(1+eta)*(1+zeta), +(1-xi)*(1+zeta), +(1-xi)*(1+eta)]], float)
    return N, dN
def solve_heat3d_h8_transient(env):
    nodes=np.asarray(env.get('nodes'),float); elems=np.asarray(env.get('elems'),int)
    k=float(env.get('k')); rho=float(env.get('rho')); cp=float(env.get('cp'))
    dt=float(env.get('dt')); t_end=float(env.get('t_end')); theta=float(env.get('theta',0.5))
    T0=env.get('T0', 0.0)
    nn=nodes.shape[0]; ndof=nn
    K=np.zeros((ndof,ndof)); M=np.zeros((ndof,ndof)); F=np.zeros(ndof)
    gp=[-1/np.sqrt(3),1/np.sqrt(3)]
    for conn in elems:
        xe=nodes[conn,:]; Ke=np.zeros((8,8)); Me=np.zeros((8,8)); Fe=np.zeros(8)
        for xi in gp:
            for eta in gp:
                for ze in gp:
                    N,dNdxi=_H8_shape(xi,eta,ze); J=dNdxi.T @ xe; detJ=np.linalg.det(J)
                    if detJ<=0: raise ValueError('Invalid H8 geometry')
                    invJ=np.linalg.inv(J); dNdx=dNdxi @ invJ.T
                    B=dNdx
                    Ke += k * (B @ B.T) * detJ
                    Me += rho*cp * (np.outer(N,N)) * detJ
                    if 'q' in env:
                        Fe += np.asarray(env['q'],float) * N * detJ
        dofs=[int(n) for n in conn]
        for i,I in enumerate(dofs):
            F[I]+=Fe[i]
            for j,J in enumerate(dofs):
                K[I,J]+=Ke[i,j]; M[I,J]+=Me[i,j]
    fixed={}
    for ent in env.get('fixT', []):
        if len(ent)==2:
            n=int(ent[0]); val=float(ent[1]); fixed[n]=val
    fixed_idx=np.array(sorted(fixed.keys()),dtype=int); free=np.setdiff1d(np.arange(ndof), fixed_idx)
    T=np.zeros(ndof) + (T0*np.ones(ndof) if np.ndim(T0)==0 else np.asarray(T0,float).reshape(ndof))
    for i in fixed_idx: T[i]=fixed[i]
    nstep=int(np.ceil(t_end/dt)); times=np.linspace(0, nstep*dt, nstep+1)
    Thist=np.zeros((nstep+1, ndof)); Thist[0,:]=T
    A = M + theta*dt*K
    for n in range(nstep):
        b = (M - (1-theta)*dt*K) @ T + dt*F
        if fixed_idx.size>0:
            b = b[free] - A[np.ix_(free, fixed_idx)] @ np.array([fixed[i] for i in fixed_idx],float)
            Tfree = np.linalg.solve(A[np.ix_(free,free)], b)
            T[free]=Tfree; T[fixed_idx]=np.array([fixed[i] for i in fixed_idx],float)
        else:
            T = np.linalg.solve(A, b)
        Thist[n+1,:]=T
    return {'T_hist':Thist, 'times':times, 'nodes':nodes}
