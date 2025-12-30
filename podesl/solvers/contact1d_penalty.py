from __future__ import annotations
import numpy as np

def solve_contact1d_penalty(env):
    L=float(env.get('L')); E=float(env.get('E')); A=float(env.get('A'))
    g=float(env.get('g',0.0)); P=float(env.get('P',0.0)); kpen=float(env.get('kpen',1e12))
    k = E*A/L
    u1_free = P/k
    if u1_free <= g:
        u1 = u1_free; R = np.array([0.0, 0.0]); U = np.array([0.0, u1])
        return {'U':U, 'contact_active': False, 'R':R, 'k':k}
    u1 = (P + kpen*g)/(k + kpen)
    R = np.array([0.0, k*u1 - P + kpen*(u1 - g)])
    U = np.array([0.0, u1])
    return {'U':U, 'contact_active': True, 'R':R, 'k':k}
