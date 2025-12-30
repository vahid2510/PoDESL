
from __future__ import annotations
from typing import Dict
import numpy as np
from ..linalg_utils import safe_solve

def _D_iso_3d(E, nu):
    lam = E*nu/((1+nu)*(1-2*nu))
    mu  = E/(2*(1+nu))
    D = np.array([
        [lam+2*mu, lam,      lam,      0,   0,   0],
        [lam,      lam+2*mu, lam,      0,   0,   0],
        [lam,      lam,      lam+2*mu, 0,   0,   0],
        [0,        0,        0,        mu,  0,   0],
        [0,        0,        0,        0,   mu,  0],
        [0,        0,        0,        0,   0,   mu],
    ], dtype=float)
    return D

def _tet4_B_matrix(xe):
    x = xe[:,0]; y = xe[:,1]; z = xe[:,2]
    M = np.array([
        [1, x[0], y[0], z[0]],
        [1, x[1], y[1], z[1]],
        [1, x[2], y[2], z[2]],
        [1, x[3], y[3], z[3]],
    ], dtype=float)
    detM = np.linalg.det(M)
    V = detM/6.0
    if V <= 0:
        raise ValueError("Invalid/negative-volume TET4")
    cof = np.linalg.inv(M).T
    gradN = cof[1:,:]  # (3,4)
    B = np.zeros((6, 12))
    for i in range(4):
        dNi_dx = gradN[0,i]
        dNi_dy = gradN[1,i]
        dNi_dz = gradN[2,i]
        B[:, 3*i:3*i+3] = np.array([
            [dNi_dx,      0.0,      0.0],
            [0.0,      dNi_dy,      0.0],
            [0.0,         0.0,   dNi_dz],
            [dNi_dy,   dNi_dx,      0.0],
            [0.0,      dNi_dz,   dNi_dy],
            [dNi_dz,      0.0,   dNi_dx],
        ])
    return B, V

def _Ke_tet4(xe, D):
    B, V = _tet4_B_matrix(xe)
    Ke = (B.T @ D @ B) * V
    return Ke, B, V

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
    ])
    dN_dxi = np.zeros((8,3))
    dN_dxi[:,0] = 0.125*np.array([
        -(1-eta)*(1-zeta),
         (1-eta)*(1-zeta),
         (1+eta)*(1-zeta),
        -(1+eta)*(1-zeta),
        -(1-eta)*(1+zeta),
         (1-eta)*(1+zeta),
         (1+eta)*(1+zeta),
        -(1+eta)*(1+zeta),
    ])
    dN_dxi[:,1] = 0.125*np.array([
        -(1-xi)*(1-zeta),
        -(1+xi)*(1-zeta),
         (1+xi)*(1-zeta),
         (1-xi)*(1-zeta),
        -(1-xi)*(1+zeta),
        -(1+xi)*(1+zeta),
         (1+xi)*(1+zeta),
         (1-xi)*(1+zeta),
    ])
    dN_dxi[:,2] = 0.125*np.array([
        -(1-xi)*(1-eta),
        -(1+xi)*(1-eta),
        -(1+xi)*(1+eta),
        -(1-xi)*(1+eta),
         (1-xi)*(1-eta),
         (1+xi)*(1-eta),
         (1+xi)*(1+eta),
         (1-xi)*(1+eta),
    ])
    return N, dN_dxi

def _hex8_B(xe, xi, eta, zeta):
    N, dN_par = _hex8_shape(xi, eta, zeta)
    J = np.zeros((3,3))
    for a in range(8):
        J[:,0] += dN_par[a,0]*xe[a,:]
        J[:,1] += dN_par[a,1]*xe[a,:]
        J[:,2] += dN_par[a,2]*xe[a,:]
    detJ = np.linalg.det(J)
    if detJ <= 0:
        raise ValueError("Invalid HEX8 mapping (detJ<=0). Check node order.")
    invJ = np.linalg.inv(J)
    dN_dx = dN_par @ invJ.T
    B = np.zeros((6, 24))
    for a in range(8):
        dNx, dNy, dNz = dN_dx[a,:]
        B[:, 3*a:3*a+3] = np.array([
            [dNx,   0.0,  0.0],
            [0.0,   dNy,  0.0],
            [0.0,   0.0,  dNz],
            [dNy,   dNx,  0.0],
            [0.0,   dNz,  dNy],
            [dNz,   0.0,  dNx],
        ])
    return B, detJ

def _Ke_hex8(xe, D):
    gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
    Ke = np.zeros((24,24))
    Blist = []
    for xi in gp:
        for eta in gp:
            for zeta in gp:
                B, detJ = _hex8_B(xe, xi, eta, zeta)
                Ke += (B.T @ D @ B) * detJ
                Blist.append(B)
    return Ke, Blist

def solve_elastic3d(env):
    nodes = np.asarray(env.get("nodes"), dtype=float)
    elems = np.asarray(env.get("elems"), dtype=int)
    etype = env.get("etype", "HEX8").upper()
    E = float(env.get("E", 210e9))
    nu = float(env.get("nu", 0.3))
    D = _D_iso_3d(E, nu)

    nn = nodes.shape[0]; ndof = 3*nn
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)
    elem_data = []
    for e, conn in enumerate(elems):
        xe = nodes[conn,:]
        if etype=="TET4":
            if conn.size!=4: raise ValueError("TET4 needs 4 nodes")
            Ke, B, V = _Ke_tet4(xe, D)
            dofs = []
            for a in conn: dofs += [3*a, 3*a+1, 3*a+2]
            for i,I in enumerate(dofs):
                for j,J in enumerate(dofs):
                    K[I,J] += Ke[i,j]
            elem_data.append({"type":"TET4","Blist":[B],"conn":conn})
        elif etype=="HEX8":
            if conn.size!=8: raise ValueError("HEX8 needs 8 nodes")
            Ke, Blist = _Ke_hex8(xe, D)
            dofs = []
            for a in conn: dofs += [3*a, 3*a+1, 3*a+2]
            for i,I in enumerate(dofs):
                for j,J in enumerate(dofs):
                    K[I,J] += Ke[i,j]
            elem_data.append({"type":"HEX8","Blist":Blist,"conn":conn})
        else:
            raise ValueError("etype must be 'TET4' or 'HEX8'")

    for node, fx, fy, fz in env.get("point_loads", []):
        F[3*node+0] += fx; F[3*node+1] += fy; F[3*node+2] += fz

    if "body" in env:
        bx, by, bz = env["body"]
        rho = float(env.get("rho", 1.0))
        if etype=="TET4":
            for d in elem_data:
                if d["type"]!="TET4": continue
                conn = d["conn"]; xe = nodes[conn,:]
                _, V = _tet4_B_matrix(xe)
                fe_node = rho * V / 4.0
                for a in conn:
                    F[3*a+0] += fe_node*bx
                    F[3*a+1] += fe_node*by
                    F[3*a+2] += fe_node*bz
        else:
            for d in elem_data:
                if d["type"]!="HEX8": continue
                conn = d["conn"]; xe = nodes[conn,:]
                gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
                Vol = 0.0
                for xi in gp:
                    for eta in gp:
                        for zeta in gp:
                            _, detJ = _hex8_B(xe, xi, eta, zeta)
                            Vol += detJ
                fe_node = rho * Vol / 8.0
                for a in conn:
                    F[3*a+0] += fe_node*bx
                    F[3*a+1] += fe_node*by
                    F[3*a+2] += fe_node*bz

    fixed_map: Dict[int, float] = {}
    for entry in env.get("fix", []):
        if not entry:
            continue
        node = int(entry[0])
        if len(entry) >= 2:
            dof = str(entry[1]).lower()
        else:
            dof = "ux"
        val = float(entry[2]) if len(entry) >= 3 else 0.0
        if dof == "ux":
            fixed_map[3*node+0] = val
        elif dof == "uy":
            fixed_map[3*node+1] = val
        elif dof == "uz":
            fixed_map[3*node+2] = val
        elif dof == "both":
            fixed_map[3*node+0] = val
            fixed_map[3*node+1] = val
            fixed_map[3*node+2] = val
        else:
            raise ValueError("fix dof must be 'ux','uy','uz'")
    fixed_idx = np.array(sorted(fixed_map.keys()), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed_idx)

    U = np.zeros(ndof)
    if free.size>0:
        rhs = F[free]
        if fixed_idx.size > 0:
            rhs -= K[np.ix_(free, fixed_idx)] @ np.array(
                [fixed_map[i] for i in fixed_idx], float
            )
        U[free] = safe_solve(K[np.ix_(free, free)], rhs)
    if fixed_idx.size > 0:
        for idx in fixed_idx:
            U[idx] = fixed_map.get(idx, 0.0)

    sigma_node = np.zeros((nn,6)); w_node = np.zeros(nn)
    for d in elem_data:
        conn = d["conn"]
        ue = np.zeros(3*len(conn))
        for i,a in enumerate(conn):
            ue[3*i:3*i+3] = U[3*a:3*a+3]
        s_avg = np.zeros(6)
        for B in d["Blist"]:
            s_avg += (B @ ue).reshape(6)
        s_avg /= len(d["Blist"])
        for a in conn:
            sigma_node[a,:] += s_avg; w_node[a]+=1.0
    w_node[w_node==0]=1.0
    sigma_node /= w_node[:,None]

    sx, sy, sz, txy, tyz, tzx = sigma_node.T
    sigma_vm = np.sqrt(0.5*((sx-sy)**2 + (sy-sz)**2 + (sz-sx)**2) + 3.0*(txy**2 + tyz**2 + tzx**2))

    return {"U": U, "nodes": nodes, "elems": elems, "etype": etype,
            "sigma": sigma_node, "sigma_vm": sigma_vm}
