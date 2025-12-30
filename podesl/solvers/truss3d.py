from __future__ import annotations
import numpy as np
from ..linalg_utils import generalized_eig, safe_solve

def _dircos(x1,y1,z1,x2,y2,z2):
    L=float(np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2))
    if L<=0: raise ValueError('Zero-length truss element')
    l=(x2-x1)/L; m=(y2-y1)/L; n=(z2-z1)/L
    return l,m,n,L
def _k_elem_global(EA,L,l,m,n):
    gg=np.array([[l*l,l*m,l*n],[m*l,m*m,m*n],[n*l,n*m,n*n]],float)
    k3=(EA/L)*gg; k=np.zeros((6,6),float)
    k[:3,:3]+=k3; k[:3,3:]-=k3; k[3:,:3]-=k3; k[3:,3:]+=k3; return k
def solve_truss3d_static(env):
    nodes=np.asarray(env.get('nodes'),float); elems=np.asarray(env.get('elems'),int)
    Ein=env.get('E'); Ain=env.get('A'); ne=elems.shape[0]
    E=np.full(ne,float(Ein)) if np.ndim(Ein)==0 else np.asarray(Ein,float).reshape(ne)
    A=np.full(ne,float(Ain)) if np.ndim(Ain)==0 else np.asarray(Ain,float).reshape(ne)
    nn=nodes.shape[0]; ndof=3*nn; K=np.zeros((ndof,ndof)); F=np.zeros(ndof)
    for ee,(n1,n2) in enumerate(elems):
        x1,y1,z1=nodes[n1]; x2,y2,z2=nodes[n2]; l,m,n,L=_dircos(x1,y1,z1,x2,y2,z2)
        kG=_k_elem_global(E[ee]*A[ee],L,l,m,n)
        dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2]
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs): K[I,J]+=kG[i,j]
    for ent in env.get('loads',[]):
        n=int(ent[0]); fx=float(ent[1]); fy=float(ent[2]); fz=float(ent[3])
        F[3*n+0]+=fx; F[3*n+1]+=fy; F[3*n+2]+=fz
    fixed={}
    for ent in env.get('fix',[]):
        if len(ent)==3:
            n=int(ent[0]); dof=str(ent[1]).lower(); val=float(ent[2])
            if dof in ('ux','both'): fixed[3*n+0]=val
            if dof in ('uy','both'): fixed[3*n+1]=val
            if dof in ('uz','both'): fixed[3*n+2]=val
        elif len(ent)==2:
            n=int(ent[0]); val=float(ent[1]); fixed[3*n+0]=val; fixed[3*n+1]=val; fixed[3*n+2]=val
    fixed_idx=np.array(sorted(fixed.keys()),dtype=int); free_idx=np.setdiff1d(np.arange(ndof),fixed_idx)
    vals=np.zeros(ndof); 
    for i in fixed_idx: vals[i]=fixed[i]
    U=np.zeros(ndof)
    if free_idx.size>0:
        rhs=F[free_idx]-K[np.ix_(free_idx,fixed_idx)] @ vals[fixed_idx]
        U[free_idx]=safe_solve(K[np.ix_(free_idx,free_idx)], rhs)
    U[fixed_idx]=vals[fixed_idx]; R=K @ U - F
    return {'U':U,'R':R,'K':K}
def solve_truss3d_modal(env):
    nodes=np.asarray(env.get('nodes'),float); elems=np.asarray(env.get('elems'),int)
    Ein=env.get('E'); Ain=env.get('A'); rho=float(env.get('rho')); ne=elems.shape[0]
    E=np.full(ne,float(Ein)) if np.ndim(Ein)==0 else np.asarray(Ein,float).reshape(ne)
    A=np.full(ne,float(Ain)) if np.ndim(Ain)==0 else np.asarray(Ain,float).reshape(ne)
    nn=nodes.shape[0]; ndof=3*nn; K=np.zeros((ndof,ndof)); M=np.zeros((ndof,ndof))
    for ee,(n1,n2) in enumerate(elems):
        x1,y1,z1=nodes[n1]; x2,y2,z2=nodes[n2]; l,m,n,L=_dircos(x1,y1,z1,x2,y2,z2)
        kG=_k_elem_global(E[ee]*A[ee],L,l,m,n)
        dofs=[3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2]
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs): K[I,J]+=kG[i,j]
        m=rho*A[ee]*L
        for nN in (n1,n2):
            M[3*nN+0,3*nN+0]+=m/6.0; M[3*nN+1,3*nN+1]+=m/6.0; M[3*nN+2,3*nN+2]+=m/6.0
    fixed=[]
    for ent in env.get('fix',[]):
        if len(ent)==3:
            n=int(ent[0]); dof=str(ent[1]).lower()
            if dof in ('ux','both'): fixed.append(3*n+0)
            if dof in ('uy','both'): fixed.append(3*n+1)
            if dof in ('uz','both'): fixed.append(3*n+2)
        elif len(ent)==2:
            n=int(ent[0]); fixed.extend([3*n,3*n+1,3*n+2])
    fixed=np.array(sorted(set(fixed)),dtype=int); free=np.setdiff1d(np.arange(ndof),fixed)
    Kf=K[np.ix_(free,free)]; Mf=M[np.ix_(free,free)]
    lam,vec = generalized_eig(Kf, Mf)
    order=np.argsort(np.real(lam)); lam=np.real(lam[order]); vec=np.real(vec[:,order])
    omegas=np.sqrt(np.clip(lam,0,None)); freq=omegas/(2*np.pi)
    nm=int(env.get('nmodes',min(6,vec.shape[1]))); modes=np.zeros((ndof,nm))
    for k in range(nm):
        full=np.zeros(ndof); full[free]=vec[:,k]; mk=full @ (M @ full)
        if mk>0: full/=np.sqrt(mk); modes[:,k]=full
    return {'freq':freq[:nm],'omega':omegas[:nm],'modes':modes,'K':K,'M':M}
