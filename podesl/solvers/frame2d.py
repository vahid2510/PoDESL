from __future__ import annotations

import numpy as np

from ..linalg_utils import generalized_eig, safe_solve
from .bc_utils import FRAME2D_DOF_MAP, dirichlet_vector
from .input_utils import merged_dirichlet, merged_point_loads

def _k_local(EA, EI, L):
    L2=L*L; L3=L2*L
    k=np.zeros((6,6),float)
    k[0,0]=k[3,3]=EA/L; k[0,3]=k[3,0]=-EA/L
    k[1,1]=k[4,4]=12*EI/L3; k[1,4]=k[4,1]=-12*EI/L3
    k[1,2]=k[2,1]=6*EI/L2; k[1,5]=k[5,1]=6*EI/L2
    k[2,2]=k[5,5]=4*EI/L;  k[2,5]=k[5,2]=2*EI/L
    k[2,4]=k[4,2]=-6*EI/L2; k[4,5]=k[5,4]=-6*EI/L2
    return k
def _T(c,s):
    R=np.array([[c,s,0],[-s,c,0],[0,0,1]],float); T=np.zeros((6,6),float); T[:3,:3]=R; T[3:,3:]=R; return T
def _geom(x1,y1,x2,y2):
    L=float(np.hypot(x2-x1,y2-y1)); 
    if L<=0: raise ValueError('Zero-length frame element')
    c=(x2-x1)/L; s=(y2-y1)/L; return c,s,L
def _element_array(value, count):
    if value is None:
        raise KeyError("Material property missing for frame analysis.")
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return np.full(count, float(arr))
    arr = arr.reshape(-1)
    if arr.size == 1:
        return np.full(count, float(arr[0]))
    if arr.size != count:
        raise ValueError(f"Property array must have length {count}, got {arr.size}.")
    return arr
def solve_frame2d_static(env):
    nodes = np.asarray(env.get("nodes"), float)
    elems = np.asarray(env.get("elems"), int)
    ne = elems.shape[0]
    nn = nodes.shape[0]
    ndof = 3 * nn

    E = _element_array(env.get("E"), ne)
    A = _element_array(env.get("A"), ne)
    I = _element_array(env.get("I"), ne)

    K = np.zeros((ndof, ndof), float)
    F = np.zeros(ndof, float)

    csL = []

    for ee, (n1, n2) in enumerate(elems):
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]
        c, s, L = _geom(x1, y1, x2, y2)
        csL.append((c, s, L))
        kL = _k_local(E[ee] * A[ee], E[ee] * I[ee], L)
        T = _T(c, s)
        kG = T.T @ kL @ T
        dofs = np.array([3 * n1, 3 * n1 + 1, 3 * n1 + 2, 3 * n2, 3 * n2 + 1, 3 * n2 + 2])
        for i, Iglob in enumerate(dofs):
            for j, Jglob in enumerate(dofs):
                K[Iglob, Jglob] += kG[i, j]

    load_entries = merged_point_loads(
        env.get("loads"),
        env.get("point_loads"),
        dof_order=("ux", "uy", "rz"),
    )
    for node, comp in load_entries:
        F[3 * node + 0] += comp[0]
        F[3 * node + 1] += comp[1]
        F[3 * node + 2] += comp[2] if len(comp) > 2 else 0.0

    bc_entries = merged_dirichlet(
        env.get("fix"),
        env.get("bcs"),
        dof_map=FRAME2D_DOF_MAP,
    )
    fixed_idx, fixed_vals = dirichlet_vector(
        bc_entries,
        ndof=ndof,
        ndof_per_node=3,
        dof_map=FRAME2D_DOF_MAP,
    )
    free_idx = np.setdiff1d(np.arange(ndof), fixed_idx)

    U = np.zeros(ndof)
    if free_idx.size > 0:
        rhs = F[free_idx]
        if fixed_idx.size:
            rhs -= K[np.ix_(free_idx, fixed_idx)] @ fixed_vals[fixed_idx]
        U[free_idx] = safe_solve(K[np.ix_(free_idx, free_idx)], rhs)
    if fixed_idx.size:
        U[fixed_idx] = fixed_vals[fixed_idx]

    R = K @ U - F

    member_forces = np.zeros((ne, 6))
    for ee, (n1, n2) in enumerate(elems):
        c, s, L = csL[ee]
        T = _T(c, s)
        kL = _k_local(E[ee] * A[ee], E[ee] * I[ee], L)
        dofs = np.array([3 * n1, 3 * n1 + 1, 3 * n1 + 2, 3 * n2, 3 * n2 + 1, 3 * n2 + 2])
        u_local = T @ U[dofs]
        member_forces[ee, :] = kL @ u_local

    disp = U.reshape(nn, 3)

    return {
        "U": U,
        "disp": disp,
        "R": R,
        "reactions": R,
        "K": K,
        "member_forces": member_forces,
        "nodes": nodes,
        "elems": elems,
    }
def solve_frame2d_modal(env):
    nodes=np.asarray(env.get('nodes'),float); elems=np.asarray(env.get('elems'),int)
    E=float(env.get('E')); A=float(env.get('A')); I=float(env.get('I')); rho=float(env.get('rho'))
    ndof=3*nodes.shape[0]; K=np.zeros((ndof,ndof),float); M=np.zeros((ndof,ndof),float)
    for (n1,n2) in elems:
        x1,y1=nodes[n1]; x2,y2=nodes[n2]; c,s,L=_geom(x1,y1,x2,y2)
        kL=_k_local(E*A,E*I,L); T=_T(c,s); kG=T.T @ kL @ T
        dofs=np.array([3*n1,3*n1+1,3*n1+2,3*n2,3*n2+1,3*n2+2])
        m=rho*A*L; J=rho*I*L
        m_local=np.diag([m/2,m/2,J/2,m/2,m/2,J/2]); mG=T.T @ m_local @ T
        for i, Iglob in enumerate(dofs):
            for j, Jglob in enumerate(dofs):
                K[Iglob,Jglob]+=kG[i,j]; M[Iglob,Jglob]+=mG[i,j]
    fixed=[]
    for ent in env.get('fix',[]):
        if len(ent)==3:
            n=int(ent[0]); dof=str(ent[1]).lower()
            if dof in ('ux','both'): fixed.append(3*n+0)
            if dof in ('uy','both'): fixed.append(3*n+1)
            if dof in ('rz','both'): fixed.append(3*n+2)
        elif len(ent)==2:
            n=int(ent[0]); fixed.extend([3*n+0,3*n+1,3*n+2])
    fixed=np.array(sorted(set(fixed)),dtype=int); free=np.setdiff1d(np.arange(ndof),fixed)
    Kf=K[np.ix_(free,free)]; Mf=M[np.ix_(free,free)]
    lam,vec = generalized_eig(Kf, Mf)
    order=np.argsort(np.real(lam)); lam=np.real(lam[order]); vec=np.real(vec[:,order])
    omegas=np.sqrt(np.clip(lam,0,None)); freq=omegas/(2*np.pi)
    nm=int(env.get('nmodes',min(6,vec.shape[1]))); modes=np.zeros((ndof,nm))
    for k in range(nm):
        full=np.zeros(ndof); full[free]=vec[:,k]; mk=full @ (M @ full); 
        if mk>0: full/=np.sqrt(mk); modes[:,k]=full
    return {'freq':freq[:nm],'omega':omegas[:nm],'modes':modes,'K':K,'M':M}
