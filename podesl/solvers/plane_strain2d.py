
from __future__ import annotations
import numpy as np
from .bc_utils import PLANAR_UV_DOF_MAP, dirichlet_vector

def _D_plane_strain(E, nu):
    c = E / ((1.0 + nu) * (1.0 - 2.0*nu))
    D = c * np.array([
        [1.0 - nu,     nu,           0.0],
        [    nu,   1.0 - nu,         0.0],
        [   0.0,       0.0, (1.0 - 2.0*nu)/2.0]
    ], dtype=float)
    return D

def _K_T3(xe, D, t):
    x1,y1 = xe[0]; x2,y2 = xe[1]; x3,y3 = xe[2]
    A = 0.5*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    if A == 0: raise ValueError("Zero-area triangle")
    b = np.array([y2-y3, y3-y1, y1-y2], dtype=float)
    c = np.array([x3-x2, x1-x3, x2-x1], dtype=float)
    B = (1.0/(2*A))*np.array([
        [b[0], 0.0, b[1], 0.0, b[2], 0.0],
        [0.0, c[0], 0.0, c[1], 0.0, c[2]],
        [c[0], b[0], c[1], b[1], c[2], b[2]]
    ], dtype=float)
    Ke = (B.T @ D @ B) * (2*A) * t
    return Ke, B, A

def _K_Q4(xe, D, t):
    gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
    Ke = np.zeros((8,8)); Blist=[]
    for xi in gp:
        for eta in gp:
            Nxi = 0.25*np.array([-(1-eta),  (1-eta),  (1+eta), -(1+eta)])
            Neta= 0.25*np.array([-(1-xi),  -(1+xi),   (1+xi),   (1-xi)])
            J = np.zeros((2,2))
            for a in range(4):
                J[0,0] += Nxi[a]*xe[a,0]; J[0,1] += Nxi[a]*xe[a,1]
                J[1,0] += Neta[a]*xe[a,0];J[1,1] += Neta[a]*xe[a,1]
            detJ = np.linalg.det(J)
            if detJ <= 0: raise ValueError("Invalid Q4 element mapping (detJ<=0)")
            invJ = np.linalg.inv(J)
            dNdx = np.zeros((4,2))
            for a in range(4):
                grad = invJ @ np.array([Nxi[a], Neta[a]])
                dNdx[a,:] = grad
            B = np.zeros((3,8))
            for a in range(4):
                B[0,2*a]   = dNdx[a,0]
                B[1,2*a+1] = dNdx[a,1]
                B[2,2*a]   = dNdx[a,1]
                B[2,2*a+1] = dNdx[a,0]
            Ke += (B.T @ D @ B) * detJ * t
            Blist.append(B)
    return Ke, Blist

def solve_plane_strain2d(env):
    nodes = np.asarray(env.get('nodes'), dtype=float)
    elems = np.asarray(env.get('elems'), dtype=int)
    etype = env.get('etype', 'Q4').upper()
    E = float(env.get('E', 210e9))
    nu = float(env.get('nu', 0.3))
    t = float(env.get('t', 1.0))

    nn = nodes.shape[0]; ndof = 2*nn
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    D = _D_plane_strain(E, nu)

    elem_data = []
    for e,conn in enumerate(elems):
        xe = nodes[conn,:]
        if etype=='T3' and conn.size!=3:
            raise ValueError('T3 needs 3 nodes per element')
        if etype=='Q4' and conn.size!=4:
            raise ValueError('Q4 needs 4 nodes per element')
        if etype=='T3':
            Ke, B, A = _K_T3(xe, D, t)
            dofs = []
            for a in conn: dofs += [2*a, 2*a+1]
            for i,I in enumerate(dofs):
                for j,J in enumerate(dofs):
                    K[I,J] += Ke[i,j]
            elem_data.append({'type':'T3','B':B,'A':A,'conn':conn})
        else:
            Ke, Blist = _K_Q4(xe, D, t)
            dofs = []
            for a in conn: dofs += [2*a, 2*a+1]
            for i,I in enumerate(dofs):
                for j,J in enumerate(dofs):
                    K[I,J] += Ke[i,j]
            elem_data.append({'type':'Q4','Blist':Blist,'conn':conn})

    for node, fx, fy in env.get('point_loads', []):
        F[2*node] += fx; F[2*node+1] += fy

    if 'body' in env:
        bx, by = env['body']
        rho = float(env.get('rho', 1.0))
        gfac = rho
        if etype=='T3':
            for d in elem_data:
                if d['type']!='T3': continue
                conn = d['conn']; A = d['A']
                fe = (A * t * gfac) * np.array([bx,by,bx,by,bx,by], dtype=float)/3.0
                dofs = []
                for a in conn: dofs += [2*a, 2*a+1]
                F[dofs] += fe
        else:
            for d in elem_data:
                if d['type']!='Q4': continue
                conn = d['conn']
                xe = nodes[conn,:]
                Aapprox = 0.5*abs(np.linalg.det(np.array([xe[1]-xe[0], xe[3]-xe[0]]))) +                           0.5*abs(np.linalg.det(np.array([xe[1]-xe[2], xe[3]-xe[2]])))
                fe_node = (Aapprox * t * gfac) / 4.0
                for a in conn:
                    F[2*a] += fe_node*bx; F[2*a+1] += fe_node*by

    for tr in env.get('traction', []):
        if isinstance(tr, dict):
            ee = int(tr['elem']); edge = int(tr['edge']); tx = float(tr.get('tx',0.0)); ty = float(tr.get('ty',0.0))
        else:
            ee, edge, tx, ty = int(tr[0]), int(tr[1]), float(tr[2]), float(tr[3])
        conn = elems[ee]
        if conn.size == 3: edge_nodes = {0:(0,1), 1:(1,2), 2:(2,0)}[edge]
        else:              edge_nodes = {0:(0,1), 1:(1,2), 2:(2,3), 3:(3,0)}[edge]
        n1 = conn[edge_nodes[0]]; n2 = conn[edge_nodes[1]]
        x1,y1 = nodes[n1]; x2,y2 = nodes[n2]
        Ledge = float(np.hypot(x2-x1, y2-y1))
        fe = (Ledge * t) * 0.5
        F[2*n1]   += tx * fe; F[2*n1+1] += ty * fe
        F[2*n2]   += tx * fe; F[2*n2+1] += ty * fe

    fixed_idx, fixed_vals = dirichlet_vector(
        env.get('fix', []),
        ndof=ndof,
        ndof_per_node=2,
        dof_map=PLANAR_UV_DOF_MAP,
    )
    free = np.setdiff1d(np.arange(ndof), fixed_idx)

    U = np.zeros(ndof)
    if free.size>0:
        rhs = F[free].copy()
        if fixed_idx.size>0:
            rhs -= K[np.ix_(free, fixed_idx)] @ fixed_vals[fixed_idx]
        U[free] = np.linalg.solve(K[np.ix_(free, free)], rhs)
    if fixed_idx.size:
        U[fixed_idx] = fixed_vals[fixed_idx]

    sigma_elem = []
    sigma_node = np.zeros((nn,3)); w_node = np.zeros(nn)
    for d in elem_data:
        conn = d['conn']
        ue = np.zeros(2*len(conn))
        for i,a in enumerate(conn):
            ue[2*i:2*i+2] = U[2*a:2*a+2]
        if d['type']=='T3':
            B = d['B']
            s = (B @ ue).reshape(3)
            sigma_elem.append(s)
            for a in conn:
                sigma_node[a,:] += s; w_node[a]+=1.0
        else:
            Blist = d['Blist']
            s_avg = np.zeros(3)
            for B in Blist:
                s_avg += (B @ ue).reshape(3)
            s_avg /= len(Blist)
            sigma_elem.append(s_avg)
            for a in conn:
                sigma_node[a,:] += s_avg; w_node[a]+=1.0
    w_node[w_node==0]=1.0
    sigma_node /= w_node[:,None]
    sigma_elem = np.vstack(sigma_elem) if sigma_elem else np.zeros((0,3))

    sx = sigma_node[:,0]; sy = sigma_node[:,1]; txy = sigma_node[:,2]
    sz = nu * (sx + sy)
    sigma_vm = np.sqrt(
        0.5*((sx-sy)**2 + (sy-sz)**2 + (sz-sx)**2) + 3.0*(txy**2)
    )

    return {'U': U, 'nodes': nodes, 'elems': elems, 'etype': etype,
            'sigma': sigma_node, 'sigma_vm': sigma_vm, 'sigma_elem': sigma_elem}
