from __future__ import annotations
import numpy as np

def _T(c,s):
    R=np.array([[c,s,0],[-s,c,0],[0,0,1]],float); T=np.zeros((6,6),float); T[:3,:3]=R; T[3:,3:]=R; return T
def _geom(x1,y1,x2,y2):
    L=float(np.hypot(x2-x1,y2-y1)); 
    if L<=0: raise ValueError('Zero-length frame element')
    c=(x2-x1)/L; s=(y2-y1)/L; return c,s,L
def _k_local(EA,EI,L):
    L2=L*L; L3=L2*L; k=np.zeros((6,6),float)
    k[0,0]=k[3,3]=EA/L; k[0,3]=k[3,0]=-EA/L
    k[1,1]=k[4,4]=12*EI/L3; k[1,4]=k[4,1]=-12*EI/L3
    k[1,2]=k[2,1]=6*EI/L2;  k[1,5]=k[5,1]=6*EI/L2
    k[2,2]=k[5,5]=4*EI/L;   k[2,5]=k[5,2]=2*EI/L
    k[2,4]=k[4,2]=-6*EI/L2; k[4,5]=k[5,4]=-6*EI/L2
    return k
def _fe_udl_local(qx,qy,L):
    fx1= qx*L/2; fx2= qx*L/2
    fy1= qy*L/2; fy2= qy*L/2
    m1 =  qy*L*L/12; m2 = -qy*L*L/12
    return np.array([fx1, fy1, m1, fx2, fy2, m2], float)
def _apply_fix(ndof,fix_list):
    fixed={}
    for ent in fix_list:
        if len(ent)==3:
            n=int(ent[0]); dof=str(ent[1]).lower(); val=float(ent[2])
            if dof in ('ux','both'): fixed[3*n+0]=val
            if dof in ('uy','both'): fixed[3*n+1]=val
            if dof in ('rz','both'): fixed[3*n+2]=val
        elif len(ent)==2:
            n=int(ent[0]); val=float(ent[1]); fixed[3*n]=val; fixed[3*n+1]=val; fixed[3*n+2]=val
    fixed_idx=np.array(sorted(fixed.keys()),dtype=int); vals=np.zeros(ndof)
    for i in fixed_idx: vals[i]=fixed[i]
    return fixed_idx, vals
def solve_frame2d_static(env):
    nodes=np.asarray(env.get('nodes'),float); elems=np.asarray(env.get('elems'),int)
    E=float(env.get('E')); A=float(env.get('A')); I=float(env.get('I'))
    ndof=3*nodes.shape[0]; K=np.zeros((ndof,ndof)); F=np.zeros(ndof)
    for ee,(n1,n2) in enumerate(elems):
        x1,y1=nodes[n1]; x2,y2=nodes[n2]; c,s,L=_geom(x1,y1,x2,y2)
        kL=_k_local(E*A,E*I,L); T=_T(c,s); kG=T.T @ kL @ T
        dofs=np.array([3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2])
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs): K[I,J]+=kG[i,j]
        for spec in env.get('edist', []):
            if int(spec[0])==ee:
                qx=float(spec[1]); qy=float(spec[2])
                feL=_fe_udl_local(qx,qy,L); feG=T.T @ feL
                for i,I in enumerate(dofs): F[I]+=feG[i]
    for ent in env.get('loads', []):
        n=int(ent[0]); fx=float(ent[1]); fy=float(ent[2]); mz=0.0
        if len(ent)>3: mz=float(ent[3])
        F[3*n+0]+=fx; F[3*n+1]+=fy; F[3*n+2]+=mz
    fixed_idx, vals=_apply_fix(ndof, env.get('fix',[])); free=np.setdiff1d(np.arange(ndof),fixed_idx)
    U=np.zeros(ndof)
    if free.size>0:
        rhs=F[free]-K[np.ix_(free,fixed_idx)] @ vals[fixed_idx]
        U[free]=np.linalg.solve(K[np.ix_(free,free)], rhs)
    U[fixed_idx]=vals[fixed_idx]; R=K @ U - F
    return {'U':U,'R':R,'K':K}
