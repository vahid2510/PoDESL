
from __future__ import annotations
import numpy as np
from ..linalg_utils import safe_solve

def _tet4_grad(xe):
    M = np.ones((4,4)); M[:,1:] = xe
    detM = np.linalg.det(M); V = detM/6.0
    if V <= 0: raise ValueError("Invalid/negative-volume TET4")
    cof = np.linalg.inv(M).T
    gradN = cof[1:,:]  # (3,4)
    return gradN, V

def _hex8_shape(xi, eta, zeta):
    N = np.array([
        0.125*(1-xi)*(1-eta)*(1-zeta),
        0.125*(1+xi)*(1-eta)*(1-zeta),
        0.125*(1+xi)*(1+eta)*(1-zeta),
        0.125*(1-xi)*(1+eta)*(1-zeta),
        0.125*(1-xi)*(1-eta)*(1+zeta),
        0.125*(1+xi)*(1-eta)*(1+zeta),
        0.125*(1+xi)*(1+eta)*(1+zeta),
        0.125*(1-xi)*(1+eta)*(1+zeta),
    ], dtype=float)
    dN = np.zeros((8,3))
    dN[:,0] = 0.125*np.array([
        -(1-eta)*(1-zeta),(1-eta)*(1-zeta),(1+eta)*(1-zeta),-(1+eta)*(1-zeta),
        -(1-eta)*(1+zeta),(1-eta)*(1+zeta),(1+eta)*(1+zeta),-(1+eta)*(1+zeta)])
    dN[:,1] = 0.125*np.array([
        -(1-xi)*(1-zeta),-(1+xi)*(1-zeta),(1+xi)*(1-zeta),(1-xi)*(1-zeta),
        -(1-xi)*(1+zeta),-(1+xi)*(1+zeta),(1+xi)*(1+zeta),(1-xi)*(1+zeta)])
    dN[:,2] = 0.125*np.array([
        -(1-xi)*(1-eta),-(1+xi)*(1-eta),-(1+xi)*(1+eta),-(1-xi)*(1+eta),
         (1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)])
    return N, dN

def _hex8_BJ(xe, xi, eta, zeta):
    N, dN_par = _hex8_shape(xi, eta, zeta)
    J = np.zeros((3,3))
    for a in range(8):
        J[:,0] += dN_par[a,0]*xe[a,:]
        J[:,1] += dN_par[a,1]*xe[a,:]
        J[:,2] += dN_par[a,2]*xe[a,:]
    detJ = np.linalg.det(J)
    if detJ <= 0: raise ValueError("Invalid HEX8 mapping (detJ<=0)")
    invJ = np.linalg.inv(J)
    dNdx = dN_par @ invJ.T
    return N, dNdx, detJ

def assemble_conduction(env):
    nodes = np.asarray(env.get("nodes"), dtype=float)
    elems = np.asarray(env.get("elems"), dtype=int)
    etype = env.get("etype","HEX8").upper()
    k = float(env.get("k", 45.0)); rho = float(env.get("rho", 7800.0)); c = float(env.get("c", 500.0))
    K = np.zeros((nodes.shape[0], nodes.shape[0])); M = np.zeros((nodes.shape[0], nodes.shape[0])); Fq = np.zeros(nodes.shape[0])
    if etype=="TET4":
        for conn in elems:
            xe = nodes[conn,:]; gradN, V = _tet4_grad(xe)
            Ke = (gradN.T @ gradN) * k * V
            Me = (rho*c*V/20.0) * (np.ones((4,4)) + np.eye(4))
            for i,I in enumerate(conn):
                for j,J in enumerate(conn):
                    K[I,J] += Ke[i,j]; M[I,J] += Me[i,j]
            q = float(env.get("q", 0.0))
            if q!=0.0: Fq[conn] += q * V / 4.0
    elif etype=="HEX8":
        gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
        for conn in elems:
            xe = nodes[conn,:]
            Ke = np.zeros((8,8)); Me = np.zeros((8,8)); fe = np.zeros(8)
            for xi in gp:
                for eta in gp:
                    for zeta in gp:
                        N, dNdx, detJ = _hex8_BJ(xe, xi, eta, zeta)
                        Ke += (dNdx @ dNdx.T) * k * detJ
                        Me += np.outer(N, N) * (rho*c) * detJ
                        fe += N * float(env.get("q", 0.0)) * detJ
            for i,I in enumerate(conn):
                for j,J in enumerate(conn):
                    K[I,J] += Ke[i,j]; M[I,J] += Me[i,j]
                Fq[I] += fe[i]
    else:
        raise ValueError("etype must be 'HEX8' or 'TET4'")
    return K, M, Fq

def apply_face_flux(env, Fq):
    nodes = np.asarray(env.get("nodes"), dtype=float); elems = np.asarray(env.get("elems"), dtype=int)
    etype = env.get("etype","HEX8").upper()
    for spec in env.get("qface", []):
        if isinstance(spec, dict):
            ee=int(spec["elem"]); face=int(spec["face"]); qn=float(spec["qn"])
        else:
            ee, face, qn = int(spec[0]), int(spec[1]), float(spec[2])
        conn = elems[ee]
        if etype=="TET4":
            faces = [(0,1,2),(0,1,3),(1,2,3),(0,2,3)]
            a,b,c = faces[face]; xa, xb, xc = nodes[conn[a]], nodes[conn[b]], nodes[conn[c]]
            area = 0.5*np.linalg.norm(np.cross(xb-xa, xc-xa))
            for nid in (a,b,c): Fq[conn[nid]] += qn*area/3.0
        else:
            faces = [(0,1,2,3),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]
            ids = faces[face]; xe = nodes[[conn[i] for i in ids]]
            gp = [-1/np.sqrt(3), 1/np.sqrt(3)]; contrib = np.zeros(4)
            for s in gp:
                for t in gp:
                    N = 0.25*np.array([(1-s)*(1-t),(1+s)*(1-t),(1+s)*(1+t),(1-s)*(1+t)])
                    dNs = 0.25*np.array([-(1-t), (1-t), (1+t), -(1+t)])
                    dNt = 0.25*np.array([-(1-s), -(1+s), (1+s), (1-s)])
                    ts = np.zeros(3); tt = np.zeros(3)
                    for i in range(4):
                        ts += dNs[i]*xe[i]; tt += dNt[i]*xe[i]
                    dA = np.linalg.norm(np.cross(ts, tt))
                    contrib += N * qn * dA
            for i, nid in enumerate(ids): Fq[conn[nid]] += contrib[i]
    return Fq

def solve_heat3d_transient(env):
    nodes = np.asarray(env.get("nodes"), dtype=float); elems = np.asarray(env.get("elems"), dtype=int)
    K, M, Fq = assemble_conduction(env); Fq = apply_face_flux(env, Fq)
    dt = float(env.get("dt", 0.1)); t_end = float(env.get("t_end", 1.0)); theta = float(env.get("theta", 0.5))
    nn = nodes.shape[0]
    fixed = {int(n):float(v) for (n,v) in env.get("fixT", [])}
    fixed_idx = np.array(sorted(fixed.keys()), dtype=int); free_idx = np.setdiff1d(np.arange(nn), fixed_idx)

    T = np.zeros(nn)
    if "T0" in env:
        T0 = env["T0"]
        if isinstance(T0, (list, tuple, np.ndarray)):
            T0 = np.asarray(T0, dtype=float); 
            if T0.size!=nn: raise ValueError("T0 size mismatch")
            T = T0.copy()
        else: T[:] = float(T0)
    for i in fixed_idx: T[i] = fixed[i]

    nsteps = int(np.round(t_end/dt)); times = np.linspace(0.0, nsteps*dt, nsteps+1)
    A_lhs = M + theta*dt*K; A_rhs = M - (1.0-theta)*dt*K

    store_every = int(env.get("store_every", 0))
    Thist = []
    if store_every>0:
        Thist.append(T.copy())

    if free_idx.size == 0:
        # All temperatures prescribed; no system solve needed.
        for _ in range(nsteps):
            if store_every>0 and (len(Thist)==0 or len(Thist)%store_every==0):
                Thist.append(T.copy())
        out = {"T": T, "times": times, "nodes": nodes, "elems": elems, "etype": env.get("etype","HEX8").upper()}
        if store_every>0:
            out["Thist"] = np.array(Thist)
        out["vtk"] = {"nodes": nodes, "elems": elems, "T": T.copy()}
        return out

    A11 = A_lhs[np.ix_(free_idx, free_idx)]
    A12 = A_lhs[np.ix_(free_idx, fixed_idx)]
    B11 = A_rhs[np.ix_(free_idx, free_idx)]
    B12 = A_rhs[np.ix_(free_idx, fixed_idx)]
    Ff = Fq[free_idx]

    for _ in range(nsteps):
        rhs = B11 @ T[free_idx] + B12 @ T[fixed_idx] + dt*Ff
        rhs -= A12 @ T[fixed_idx]
        T[free_idx] = safe_solve(A11, rhs)
        for i in fixed_idx: T[i] = fixed[i]
        if store_every>0 and (len(Thist)==0 or len(Thist)%store_every==0):
            Thist.append(T.copy())

    out = {"T": T, "times": times, "nodes": nodes, "elems": elems, "etype": env.get("etype","HEX8").upper()}
    if store_every>0: out["Thist"] = np.array(Thist)
    out["vtk"] = {"nodes": nodes, "elems": elems, "T": T.copy()}
    return out
