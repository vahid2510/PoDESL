from __future__ import annotations
import numpy as np
from ..linalg_utils import generalized_eig

def _geom(x1,y1,x2,y2):
    L=float(np.hypot(x2-x1,y2-y1))
    if L<=0: raise ValueError('Zero-length element')
    c=(x2-x1)/L; s=(y2-y1)/L
    return c,s,L
def _T(c,s):
    R=np.array([[c,s,0],[-s,c,0],[0,0,1]],float)
    T=np.zeros((6,6),float); T[:3,:3]=R; T[3:,3:]=R
    return T
def _k_local(EA, EI, L):
    L2=L*L; L3=L2*L
    k=np.zeros((6,6),float)
    k[0,0]=k[3,3]=EA/L; k[0,3]=k[3,0]=-EA/L
    k[1,1]=k[4,4]=12*EI/L3; k[1,4]=k[4,1]=-12*EI/L3
    k[1,2]=k[2,1]=6*EI/L2;  k[1,5]=k[5,1]=6*EI/L2
    k[2,2]=k[5,5]=4*EI/L;   k[2,5]=k[5,2]=2*EI/L
    k[2,4]=k[4,2]=-6*EI/L2; k[4,5]=k[5,4]=-6*EI/L2
    return k
def _m_local_consistent(rho, A, I, L):
    m_bar = rho*A*L/6.0 * np.array([[2,1],[1,2]], float)
    L2=L*L
    c = rho*A*L/420.0
    m_b = c * np.array([[156,   22*L,  54,   -13*L],
                        [22*L,  4*L2,  13*L, -3*L2],
                        [54,    13*L, 156,   -22*L],
                        [-13*L, -3*L2, -22*L, 4*L2]], float)
    M = np.zeros((6,6), float)
    M[0,0]+=m_bar[0,0]; M[0,3]+=m_bar[0,1]
    M[3,0]+=m_bar[1,0]; M[3,3]+=m_bar[1,1]
    idx=[1,2,4,5]
    for iI,I in enumerate(idx):
        for jJ,J in enumerate(idx):
            M[I,J]+=m_b[iI,jJ]
    return M
def solve_frame2d_modal_consistent(env):
    nodes=np.asarray(env.get('nodes'), float)
    elems=np.asarray(env.get('elems'), int)
    E=float(env.get('E')); A=float(env.get('A')); I=float(env.get('I')); rho=float(env.get('rho'))
    ndof=3*nodes.shape[0]
    K=np.zeros((ndof,ndof)); M=np.zeros((ndof,ndof))
    for (n1,n2) in elems:
        x1,y1=nodes[n1]; x2,y2=nodes[n2]
        c,s,L=_geom(x1,y1,x2,y2)
        kL=_k_local(E*A, E*I, L)
        mL=_m_local_consistent(rho, A, I, L)
        T=_T(c,s); kG=T.T @ kL @ T; mG=T.T @ mL @ T
        dofs=np.array([3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2])
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J]+=kG[i,j]; M[I,J]+=mG[i,j]
    fixed=[]
    for ent in env.get('fix', []):
        if len(ent)==3:
            n=int(ent[0]); dof=str(ent[1]).lower()
            if dof in ('ux','both'): fixed.append(3*n+0)
            if dof in ('uy','both'): fixed.append(3*n+1)
            if dof in ('rz','both'): fixed.append(3*n+2)
        elif len(ent)==2:
            n=int(ent[0]); fixed.extend([3*n,3*n+1,3*n+2])
    fixed=np.array(sorted(set(fixed)),dtype=int)
    free=np.setdiff1d(np.arange(ndof), fixed)
    Kf=K[np.ix_(free,free)]; Mf=M[np.ix_(free,free)]
    lam,vec = generalized_eig(Kf, Mf)
    order=np.argsort(np.real(lam)); lam=np.real(lam[order]); vec=np.real(vec[:,order])
    omegas=np.sqrt(np.clip(lam,0,None)); freq=omegas/(2*np.pi)
    nm=int(env.get('nmodes',min(6,vec.shape[1]))); modes=np.zeros((ndof,nm))
    for k in range(nm):
        full=np.zeros(ndof); full[free]=vec[:,k]
        mk=full @ (M @ full)
        if mk>0: full/=np.sqrt(mk)
        modes[:,k]=full
    return {'freq':freq[:nm], 'omega':omegas[:nm], 'modes':modes, 'K':K, 'M':M}
