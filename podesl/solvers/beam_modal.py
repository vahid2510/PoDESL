from __future__ import annotations
import numpy as np
from ..linalg_utils import generalized_eig
def _ke_beam(EI, L):
    L2=L*L; L3=L2*L
    return (EI/L3)*np.array([[12,6*L,-12,6*L],[6*L,4*L2,-6*L,2*L2],[-12,-6*L,12,-6*L],[6*L,2*L2,-6*L,4*L2]],float)
def _me_beam(rhoA, L):
    L2=L*L
    return (rhoA*L/420.0)*np.array([[156,22*L,54,-13*L],[22*L,4*L2,13*L,-3*L2],[54,13*L,156,-22*L],[-13*L,-3*L2,-22*L,4*L2]],float)
def _apply_bc_eig(K,M,fixed):
    fixed=np.array(sorted(set(int(i) for i in fixed)),int); ndof=K.shape[0]; free=np.setdiff1d(np.arange(ndof),fixed); 
    return K[np.ix_(free,free)], M[np.ix_(free,free)], free, fixed
def solve_beam_modal(env):
    L=float(env.get("L")); nel=int(env.get("nel")); E=float(env.get("E")); I=float(env.get("I"))
    rho=float(env.get("rho")); A=float(env.get("A")); nmodes=int(env.get("nmodes",5))
    left=str(env.get("left","clamped")).lower(); right=str(env.get("right","free")).lower()
    nnode=nel+1; x=np.linspace(0.0,L,nnode); ndof=2*nnode; K=np.zeros((ndof,ndof)); M=np.zeros((ndof,ndof))
    for e in range(nel):
        Le=x[e+1]-x[e]; ke=_ke_beam(E*I,Le); me=_me_beam(rho*A,Le); dofs=np.array([2*e,2*e+1,2*(e+1),2*(e+1)+1])
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J]+=ke[i,j]; M[I,J]+=me[i,j]
    fixed=[]
    def fix_node(n,bc):
        if bc=="clamped": fixed.extend([2*n,2*n+1])
        elif bc=="pinned": fixed.append(2*n)
        elif bc=="guided": fixed.append(2*n+1)
        elif bc=="free": pass
        else: raise ValueError("Unknown BC: "+str(bc))
    fix_node(0,left); fix_node(nnode-1,right)
    Kf,Mf,free,_=_apply_bc_eig(K,M,fixed)
    lam,vec = generalized_eig(Kf, Mf)
    order=np.argsort(np.real(lam)); lam=np.real(lam[order]); vec=np.real(vec[:,order])
    omegas=np.sqrt(np.clip(lam,0,None)); freq=omegas/(2*np.pi)
    nm=min(nmodes,vec.shape[1]); modes=np.zeros((ndof,nm))
    for k in range(nm):
        full=np.zeros(ndof); full[free]=vec[:,k]; mk=full @ (M @ full); 
        if mk>0: full/=np.sqrt(mk); modes[:,k]=full
    return {"freq":freq[:nm],"omega":omegas[:nm],"modes":modes,"x":x,"left":left,"right":right}
