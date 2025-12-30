from __future__ import annotations
import numpy as np
from ..linalg_utils import generalized_eig

# 1D Euler-Bernoulli beam modal analysis, 2 DOF per node [w, theta]
def _beam_km(EI, rhoA, L):
    L2 = L*L; L3 = L2*L
    k = (EI/L3)*np.array([
        [12,      6*L,    -12,     6*L],
        [6*L,  4*L2,    -6*L,  2*L2],
        [-12,    -6*L,     12,    -6*L],
        [6*L,  2*L2,    -6*L,  4*L2]], float)
    m = rhoA*L
    M = (m/420.0)*np.array([
        [156,      22*L,    54,     -13*L],
        [22*L,   4*L2,    13*L,    -3*L2],
        [54,      13*L,   156,     -22*L],
        [-13*L, -3*L2,  -22*L,    4*L2]], float)
    return k, M

def solve_beam1d_modal(env):
    x = np.asarray(env.get('nodes'), float).reshape(-1)
    elems = np.asarray(env.get('elems'), int)
    E = float(env.get('E')); I = float(env.get('I'))
    rho = float(env.get('rho')); A = float(env.get('A'))
    nmodes = int(env.get('nmodes', 6))

    nn = x.shape[0]; ndof = 2*nn
    K = np.zeros((ndof, ndof)); M = np.zeros((ndof, ndof))

    for e,(n1,n2) in enumerate(elems):
        x1,x2 = x[n1], x[n2]
        L = float(abs(x2-x1))
        if L <= 0:
            raise ValueError('Zero-length beam element')
        ke, me = _beam_km(E*I, rho*A, L)
        dofs = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
        for i,I in enumerate(dofs):
            for j,J in enumerate(dofs):
                K[I,J] += ke[i,j]
                M[I,J] += me[i,j]

    fixed = {}
    for ent in env.get('fix', []):
        if len(ent) == 3:
            n = int(ent[0]); dof = str(ent[1]).lower(); val = float(ent[2])
            if dof in ('w','both'):
                fixed[2*n] = val
            if dof in ('theta','both'):
                fixed[2*n+1] = val
        elif len(ent) == 2:
            n = int(ent[0]); val = float(ent[1])
            fixed[2*n] = val; fixed[2*n+1] = val
    fixed_idx = np.array(sorted(fixed.keys()), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed_idx)

    Kf = K[np.ix_(free, free)]
    Mf = M[np.ix_(free, free)]

    w, vec = generalized_eig(Kf, Mf)
    omegas = np.sqrt(np.clip(w, 0.0, None))
    freqs = omegas/(2.0*np.pi)

    nm = min(nmodes, vec.shape[1])
    modes = np.zeros((ndof, nm))
    for k in range(nm):
        full = np.zeros(ndof)
        full[free] = vec[:,k]
        mk = full @ (M @ full)
        if mk > 0:
            full /= np.sqrt(mk)
        modes[:,k] = full

    return {'freq': freqs[:nm], 'omega': omegas[:nm], 'modes': modes, 'K': K, 'M': M}
