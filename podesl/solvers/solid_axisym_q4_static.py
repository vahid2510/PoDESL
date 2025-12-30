from __future__ annotations
import numpy as np

from ..linalg_utils import safe_solve
def _Q4_shape(xi, eta):
    N = 0.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)],float)
    dN = 0.25*np.array([[-(1-eta),-(1-xi)],[+(1-eta),-(1+xi)],[+(1+eta),+(1+xi)],[-(1+eta),+(1-xi)]],float)
    return N, dN
def solve_solid_axisym_q4_static(env):
    nodes=np.asarray(env.get('nodes'),float); elems=np.asarray(env.get('elems'),int)
    E=float(env.get('E')); nu=float(env.get('nu')); c=E/((1+nu)*(1-2*nu))
    D=c*np.array([[1-nu,nu,nu,0],[nu,1-nu,nu,0],[nu,nu,1-nu,0],[0,0,0,0.5-nu]],float)
    nn=nodes.shape[0]; ndof=2*nn; K=np.zeros((ndof,ndof)); F=np.zeros(ndof)
    gp=[-1/np.sqrt(3),1/np.sqrt(3)]
    for conn in elems:
        xe=nodes[conn,:]; Ke=np.zeros((8,8)); fe=np.zeros(8)
        for xi in gp:
            for eta in gp:
                N,dNdxi=_Q4_shape(xi,eta); J=dNdxi.T @ xe; detJ=np.linalg.det(J)
                if detJ<=0: raise ValueError('Invalid Q4 geometry'); invJ=np.linalg.inv(J); dNdx=dNdxi @ invJ.T
                r_gp=float(np.dot(N, xe[:,0])); 
                if r_gp<=0: raise ValueError('Non-positive radius at GP')
                B=np.zeros((4,8))
                for a in range(4):
                    B[0,2*a+0]=dNdx[a,0]; B[1,2*a+1]=dNdx[a,1]; B[2,2*a+0]=N[a]/r_gp; B[3,2*a+0]=dNdx[a,1]; B[3,2*a+1]=dNdx[a,0]
                w=2*np.pi*r_gp*detJ
                Ke+=B.T @ D @ B * w
        dofs=[]; 
        for a in range(4): n=int(conn[a]); dofs += [2*n,2*n+1]
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs): K[I,J]+=Ke[i,j]
    for ent in env.get('loads',[]):
        n=int(ent[0]); fr=float(ent[1]); fz=float(ent[2]); F[2*n+0]+=fr; F[2*n+1]+=fz
    fixed={}
    for ent in env.get('fix',[]):
        if len(ent)==3:
            n=int(ent[0]); dof=str(ent[1]).lower(); val=float(ent[2])
            if dof in ('ur','both'): fixed[2*n+0]=val
            if dof in ('uz','both'): fixed[2*n+1]=val
        elif len(ent)==2:
            n=int(ent[0]); val=float(ent[1]); fixed[2*n]=val; fixed[2*n+1]=val
    fixed_idx=np.array(sorted(fixed.keys()),dtype=int); free_idx=np.setdiff1d(np.arange(ndof),fixed_idx)
    vals=np.zeros(ndof); 
    for i in fixed_idx: vals[i]=fixed[i]
    U=np.zeros(ndof)
    if free_idx.size>0:
        rhs=F[free_idx]-K[np.ix_(free_idx,fixed_idx)] @ vals[fixed_idx]
        U[free_idx]=safe_solve(K[np.ix_(free_idx,free_idx)], rhs)
    U[fixed_idx]=vals[fixed_idx]; R=K @ U - F
    return {'U':U,'R':R,'K':K}
