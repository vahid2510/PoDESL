from __future__ import annotations
import numpy as np
from ..linalg_utils import generalized_eig
from .plane2d_q4_static import solve_plane2d_q4_static

def _q4_area(nodes, conn):
    x1,y1 = nodes[conn[0]]
    x2,y2 = nodes[conn[1]]
    x3,y3 = nodes[conn[2]]
    x4,y4 = nodes[conn[3]]
    a1 = 0.5*abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    a2 = 0.5*abs((x4-x1)*(y3-y1) - (x3-x1)*(y4-y1))
    return a1 + a2

def solve_plane2d_q4_modal(env):
    nodes = np.asarray(env.get('nodes'), float)
    elems = np.asarray(env.get('elems'), int)
    rho = float(env.get('rho'))
    t = float(env.get('t', 1.0))
    nmodes = int(env.get('nmodes', 6))

    k_env = dict(env)
    res = solve_plane2d_q4_static(k_env)
    K = res['K']

    nn = nodes.shape[0]
    ndof = 2*nn
    M = np.zeros((ndof, ndof))

    for conn in elems:
        Ael = _q4_area(nodes, conn)
        m_el = rho * t * Ael
        m_node = m_el / 4.0
        for n in conn:
            M[2*n,2*n]     += m_node
            M[2*n+1,2*n+1] += m_node

    fixed = {}
    for ent in env.get('fix', []):
        if len(ent) == 3:
            n = int(ent[0]); dof = str(ent[1]).lower(); val = float(ent[2])
            if dof in ('ux','both'):
                fixed[2*n] = val
            if dof in ('uy','both'):
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
