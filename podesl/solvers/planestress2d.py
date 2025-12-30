
from __future__ import annotations
import numpy as np
from .bc_utils import PLANAR_UV_DOF_MAP, dirichlet_vector

def _D_pstress(E, nu):
    c = E/(1.0 - nu**2)
    return c * np.array([[1, nu, 0],
                         [nu, 1, 0],
                         [0,  0, (1-nu)/2]])

def _q4_shape(xi, eta):
    N = 0.25*np.array([(1-xi)*(1-eta),
                       (1+xi)*(1-eta),
                       (1+xi)*(1+eta),
                       (1-xi)*(1+eta)])
    dN_dxi = 0.25*np.array([[-(1-eta), -(1-xi)],
                            [ +(1-eta), -(1+xi)],
                            [ +(1+eta), +(1+xi)],
                            [ -(1+eta), +(1-xi)]])
    return N, dN_dxi

def _t3_B_A(xe):
    x1,y1 = xe[0]; x2,y2 = xe[1]; x3,y3 = xe[2]
    A = 0.5*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    if A <= 0: raise ValueError("Invalid/zero-area T3")
    b = np.array([y2-y3, y3-y1, y1-y2], dtype=float)
    c = np.array([x3-x2, x1-x3, x2-x1], dtype=float)
    B = (1.0/(2*A))*np.array([[b[0], 0, b[1], 0, b[2], 0],
                               [0, c[0], 0, c[1], 0, c[2]],
                               [c[0], b[0], c[1], b[1], c[2], b[2]]])
    return B, A

def _edge_pairs(etype):
    if etype=="Q4":
        return {0:(0,1), 1:(1,2), 2:(2,3), 3:(3,0)}
    else:
        return {0:(0,1), 1:(1,2), 2:(2,0)}

def solve_planestress2d(env):
    nodes = np.asarray(env.get("nodes"), dtype=float)
    elems = np.asarray(env.get("elems"), dtype=int)
    etype = env.get("etype","Q4").upper()
    E = float(env.get("E")); nu = float(env.get("nu")); t = float(env.get("t",1.0))
    D = _D_pstress(E, nu)

    nn = nodes.shape[0]; ndof = 2*nn
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    # assembly
    if etype=="Q4":
        gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
        for conn in elems:
            xe = nodes[conn,:]
            Ke = np.zeros((8,8)); fe = np.zeros(8)
            for xi in gp:
                for eta in gp:
                    N, dN_dxi = _q4_shape(xi, eta)
                    J = np.zeros((2,2))
                    for a in range(4):
                        J[0,0] += dN_dxi[a,0]*xe[a,0]; J[0,1] += dN_dxi[a,0]*xe[a,1]
                        J[1,0] += dN_dxi[a,1]*xe[a,0]; J[1,1] += dN_dxi[a,1]*xe[a,1]
                    detJ = np.linalg.det(J)
                    if detJ <= 0: raise ValueError("Invalid Q4 mapping (detJ<=0)")
                    invJ = np.linalg.inv(J)
                    dN_dx = dN_dxi @ invJ.T  # (4,2)
                    B = np.zeros((3,8))
                    for a in range(4):
                        B[0,2*a+0] = dN_dx[a,0]
                        B[1,2*a+1] = dN_dx[a,1]
                        B[2,2*a+0] = dN_dx[a,1]
                        B[2,2*a+1] = dN_dx[a,0]
                    Ke += (B.T @ D @ B) * t * detJ
                    if "b" in env:
                        bx, by = env["b"]
                        fe += np.kron(N, np.array([bx, by])) * t * detJ
            dofs = np.array([[2*i,2*i+1] for i in conn]).ravel()
            for i,I in enumerate(dofs):
                for j,J in enumerate(dofs):
                    K[I,J] += Ke[i,j]
            F[dofs] += fe
    elif etype=="T3":
        for conn in elems:
            xe = nodes[conn,:]
            B, A = _t3_B_A(xe)
            Ke = (B.T @ D @ B) * t * A
            fe = np.zeros(6)
            if "b" in env:
                bx, by = env["b"]
                fe = (t*A/3.0) * np.array([bx,by,bx,by,bx,by])
            dofs = np.array([[2*i,2*i+1] for i in conn]).ravel()
            for i,I in enumerate(dofs):
                for j,J in enumerate(dofs):
                    K[I,J] += Ke[i,j]
            F[dofs] += fe
    else:
        raise ValueError("etype must be Q4 or T3")

    # edge tractions
    for spec in env.get("tedge", []):
        if isinstance(spec, dict):
            ee=int(spec["elem"]); edge=int(spec["edge"]); tx=float(spec.get("tx",0.0)); ty=float(spec.get("ty",0.0))
        else:
            ee, edge, tx, ty = int(spec[0]), int(spec[1]), float(spec[2]), float(spec[3])
        conn = elems[ee]
        a,b = _edge_pairs(etype)[edge]
        n1 = conn[a]; n2 = conn[b]
        x1,y1 = nodes[n1]; x2,y2 = nodes[n2]
        L = float(np.hypot(x2-x1, y2-y1))
        F[2*n1:2*n1+2] += t * L * 0.5 * np.array([tx,ty])
        F[2*n2:2*n2+2] += t * L * 0.5 * np.array([tx,ty])

    fixed_idx, fixed_vals = dirichlet_vector(
        env.get("fix", []),
        ndof=ndof,
        ndof_per_node=2,
        dof_map=PLANAR_UV_DOF_MAP,
    )
    free_idx = np.setdiff1d(np.arange(ndof), fixed_idx)

    U = np.zeros(ndof)
    if free_idx.size>0:
        rhs = F[free_idx].copy()
        if fixed_idx.size>0:
            rhs -= K[np.ix_(free_idx, fixed_idx)] @ fixed_vals[fixed_idx]
        U[free_idx] = np.linalg.solve(K[np.ix_(free_idx, free_idx)], rhs)
    if fixed_idx.size:
        U[fixed_idx] = fixed_vals[fixed_idx]

    # element stresses
    sigma_elem = []
    sigma_node = np.zeros((nn,3)); counts = np.zeros(nn)
    if etype=="T3":
        for conn in elems:
            xe = nodes[conn,:]; B, A = _t3_B_A(xe)
            dofs = np.array([[2*i,2*i+1] for i in conn]).ravel()
            eps = B @ U[dofs]
            sig = _D_pstress(E, nu) @ eps
            sigma_elem.append(sig)
            for a in conn:
                sigma_node[a,:] += sig
                counts[a] += 1.0
    else:
        gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
        for conn in elems:
            xe = nodes[conn,:]
            sigma_sum = np.zeros(3); wsum = 0.0
            for xi in gp:
                for eta in gp:
                    N, dN_dxi = _q4_shape(xi, eta)
                    J = np.zeros((2,2))
                    for a in range(4):
                        J[0,0] += dN_dxi[a,0]*xe[a,0]; J[0,1] += dN_dxi[a,0]*xe[a,1]
                        J[1,0] += dN_dxi[a,1]*xe[a,0]; J[1,1] += dN_dxi[a,1]*xe[a,1]
                    detJ = np.linalg.det(J)
                    invJ = np.linalg.inv(J)
                    dN_dx = dN_dxi @ invJ.T
                    B = np.zeros((3,8))
                    for a in range(4):
                        B[0,2*a+0] = dN_dx[a,0]
                        B[1,2*a+1] = dN_dx[a,1]
                    B[2,2*a+0] = dN_dx[a,1]
                    B[2,2*a+1] = dN_dx[a,0]
                dofs = np.array([[2*i,2*i+1] for i in conn]).ravel()
                eps = B @ U[dofs]
                sigma_sum += (_D_pstress(E, nu) @ eps) * detJ
                wsum += detJ
            s_avg = sigma_sum/wsum
            sigma_elem.append(s_avg)
            for a in conn:
                sigma_node[a,:] += s_avg
                counts[a] += 1.0
    counts[counts==0] = 1.0
    sigma_node /= counts[:,None]
    sigma_elem_arr = np.vstack(sigma_elem) if sigma_elem else np.zeros((0,3))
    sx = sigma_node[:,0]; sy = sigma_node[:,1]; txy = sigma_node[:,2]
    sigma_vm = np.sqrt(sx*sx + sy*sy - sx*sy + 3.0*(txy**2))
    return {
        "U": U,
        "sigma": sigma_node,
        "sigma_elem": sigma_elem_arr,
        "sigma_vm": sigma_vm,
        "nodes": nodes,
        "elems": elems,
        "etype": etype,
    }
