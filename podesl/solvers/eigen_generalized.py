from __future__ import annotations
import numpy as np
from ..linalg_utils import generalized_eig
def solve_eigen_generalized(env):
    K=np.array(env.get("K"),float); M=np.array(env.get("M"),float); nmodes=int(env.get("nmodes",min(K.shape[0],6)))
    lam,vec = generalized_eig(K, M)
    order=np.argsort(np.real(lam)); lam=np.real(lam[order]); vec=np.real(vec[:,order])
    omegas=np.sqrt(np.clip(lam,0,None)); freq=omegas/(2*np.pi)
    modes=vec[:,:nmodes].copy()
    for k in range(modes.shape[1]):
        mk=modes[:,k] @ (M @ modes[:,k])
        if mk>0: modes[:,k]/=np.sqrt(mk)
    return {"lambda":lam[:nmodes],"omega":omegas[:nmodes],"freq":freq[:nmodes],"modes":modes}
