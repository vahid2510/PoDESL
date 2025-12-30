from __future__ import annotations
import numpy as np
from ..linalg_utils import generalized_eig

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
def _kg_local(N,L):
    L2=L*L; c1=N/(30.0*L)
    return c1*np.array([[0,0,0,0,0,0],
                        [0,36,3*L,0,-36,3*L],
                        [0,3*L,4*L2,0,-3*L,-1*L2],
                        [0,0,0,0,0,0],
                        [0,-36,-3*L,0,36,-3*L],
                        [0,3*L,-1*L2,0,-3*L,4*L2]],float)
def _apply_fix(ndof,fix_list):
    fixed=[]
    for ent in fix_list:
        if len(ent)==3:
            n=int(ent[0]); dof=str(ent[1]).lower()
            if dof in ('ux','both'): fixed.append(3*n+0)
            if dof in ('uy','both'): fixed.append(3*n+1)
            if dof in ('rz','both'): fixed.append(3*n+2)
        elif len(ent)==2:
            n=int(ent[0]); fixed.extend([3*n,3*n+1,3*n+2])
    return np.array(sorted(set(fixed)),dtype=int)
def solve_frame2d_buckling(env):
    nodes=np.asarray(env.get('nodes'),float); elems=np.asarray(env.get('elems'),int)
    E=float(env.get('E')); A=float(env.get('A')); I=float(env.get('I'))
    Nin=env.get('Nref',0.0); ne=elems.shape[0]
    N=np.full(ne,float(Nin)) if np.ndim(Nin)==0 else np.asarray(Nin,float).reshape(ne)
    ndof=3*nodes.shape[0]; K=np.zeros((ndof,ndof)); Kg=np.zeros((ndof,ndof))
    for ee,(n1,n2) in enumerate(elems):
        x1,y1=nodes[n1]; x2,y2=nodes[n2]; c,s,L=_geom(x1,y1,x2,y2)
        kL=_k_local(E*A,E*I,L); kgL=_kg_local(N[ee],L); T=_T(c,s)
        kG=T.T @ kL @ T; kgG=T.T @ kgL @ T
        dofs=np.array([3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2])
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J]+=kG[i,j]; Kg[I,J]+=kgG[i,j]
    fixed=_apply_fix(ndof, env.get('fix',[])); free=np.setdiff1d(np.arange(ndof),fixed)
    Kf=K[np.ix_(free,free)]; Kgf=Kg[np.ix_(free,free)]
    lam,vec = generalized_eig(Kf, -Kgf)
    order=np.argsort(np.real(lam)); lam=np.real(lam[order]); vec=np.real(vec[:,order])
    lam_cr=lam; nm=int(env.get('nmodes',min(6,vec.shape[1]))); modes=np.zeros((ndof,nm))
    for k in range(nm):
        full=np.zeros(ndof); full[free]=vec[:,k]; modes[:,k]=full/np.linalg.norm(full)
    return {'lambda_cr':lam_cr[:nm],'modes':modes,'K':K,'Kg':Kg}
