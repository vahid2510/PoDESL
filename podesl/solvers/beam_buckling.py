from __future__ import annotations
import numpy as np
def _ke_beam(EI,L):
    L2=L*L; L3=L2*L
    return (EI/L3)*np.array([[12,6*L,-12,6*L],[6*L,4*L2,-6*L,2*L2],[-12,-6*L,12,-6*L],[6*L,2*L2,-6*L,4*L2]],float)
def _kg_beam(N,L):
    L2=L*L
    return (N/(30.0*L))*np.array([[36,3*L,-36,3*L],[3*L,4*L2,-3*L,-1*L2],[-36,-3*L,36,-3*L],[3*L,-1*L, -3*L,4*L2]],float)
def _apply_bc(K,fixed):
    fixed=np.array(sorted(set(int(i) for i in fixed)),int); ndof=K.shape[0]; free=np.setdiff1d(np.arange(ndof),fixed)
    return K[np.ix_(free,free)], free, fixed
def solve_beam_buckling(env):
    L=float(env.get('L')); nel=int(env.get('nel')); E=float(env.get('E')); I=float(env.get('I')); Nref=float(env.get('Nref',1.0))
    left=str(env.get('left','pinned')).lower(); right=str(env.get('right','pinned')).lower(); nmodes=int(env.get('nmodes',4))
    nnode=nel+1; x=np.linspace(0.0,L,nnode); ndof=2*nnode; K=np.zeros((ndof,ndof)); Kg=np.zeros((ndof,ndof))
    for e in range(nel):
        Le=x[e+1]-x[e]; ke=_ke_beam(E*I,Le); kg=_kg_beam(Nref,Le); dofs=np.array([2*e,2*e+1,2*(e+1),2*(e+1)+1])
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J]+=ke[i,j]; Kg[I,J]+=kg[i,j]
    fixed=[]
    def fix_node(n,bc):
        if bc=='clamped': fixed.extend([2*n,2*n+1])
        elif bc=='pinned': fixed.append(2*n)
        elif bc=='guided': fixed.append(2*n+1)
        elif bc=='free': pass
        else: raise ValueError('Unknown BC: '+str(bc))
    fix_node(0,left); fix_node(nnode-1,right)
    Kf,free,_=_apply_bc(K,fixed); Kg_f,_,_=_apply_bc(Kg,fixed)
    try:
        lam,vec=np.linalg.eig(np.linalg.solve(Kg_f,Kf))
    except np.linalg.LinAlgError:
        lam,vec=np.linalg.eig(np.linalg.solve(Kg_f+1e-12*np.eye(Kg_f.shape[0]),Kf))
    lam=np.real(lam); order=np.argsort(lam); lam=lam[order]; vec=np.real(vec[:,order])
    Pcr=lam*Nref; nm=min(nmodes,vec.shape[1]); modes=np.zeros((ndof,nm))
    for k in range(nm):
        full=np.zeros(ndof); full[free]=vec[:,k]; 
        s=np.max(np.abs(full)); modes[:,k]=full/s if s>0 else full
    return {'lambda':lam[:nm],'Pcr':Pcr[:nm],'modes':modes,'x':x,'left':left,'right':right}
