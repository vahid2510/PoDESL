
from __future__ import annotations
import numpy as np
from .bc_utils import PLANAR_UV_DOF_MAP, dirichlet_vector

def _q4_shape(xi, eta):
    N = 0.25*np.array([(1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta)], float)
    dN = 0.25*np.array([[-(1-eta),-(1-xi)],[+(1-eta),-(1+xi)],[+(1+eta),+(1+xi)],[-(1+eta),+(1-xi)]], float)
    return N, dN

def _assemble_q4_lin(nodes, elems, D, env):
    nodes = np.asarray(nodes, float); elems = np.asarray(elems, int)
    nn = nodes.shape[0]; ndof = 2*nn
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
    for conn in elems:
        xe = nodes[conn,:]   # (4,2)
        Ke = np.zeros((8,8)); fe = np.zeros(8)
        for xi in gp:
            for eta in gp:
                N, dNdxi = _q4_shape(xi, eta)
                J = dNdxi.T @ xe; detJ = np.linalg.det(J)
                if detJ <= 0: raise ValueError('Invalid Q4 geometry')
                invJ = np.linalg.inv(J); dNdx = dNdxi @ invJ.T  # (4,2)
                B = np.zeros((3,8))
                for a in range(4):
                    B[0,2*a+0] = dNdx[a,0]        # exx
                    B[1,2*a+1] = dNdx[a,1]        # eyy
                    B[2,2*a+0] = dNdx[a,1]        # gxy
                    B[2,2*a+1] = dNdx[a,0]
                w = detJ
                Ke += B.T @ D @ B * w
                if 'bf' in env:
                    bx, by = np.asarray(env['bf'], float).tolist()
                    for a in range(4):
                        fe[2*a+0] += N[a]*bx*w
                        fe[2*a+1] += N[a]*by*w
        # edge tractions: hedge=[[elem, edge, tx, ty], ...] edges: 0-1,1-2,2-3,3-0
        for spec in env.get('hedge', []):
            if int(spec[0]) != int(conn.tolist().__hash__()): 
                pass
        # assemble to global
        dofs = []
        for a in range(4):
            n = int(conn[a]); dofs += [2*n, 2*n+1]
        for i,I in enumerate(dofs):
            F[I] += fe[i]
            for j,J in enumerate(dofs):
                K[I,J] += Ke[i,j]
    # edge tractions once more, with explicit loop over specs
    for spec in env.get('hedge', []):
        ee=int(spec[0]); edge=int(spec[1]); tx=float(spec[2]); ty=float(spec[3])
        conn = elems[ee]
        a,b = [(0,1),(1,2),(2,3),(3,0)][edge]
        n1,n2 = conn[a], conn[b]
        x1,y1 = nodes[n1]; x2,y2 = nodes[n2]
        L = float(np.hypot(x2-x1,y2-y1))
        # consistent 2-node edge load
        Fe = (L/6.0)*np.array([2*tx,2*ty, tx,ty])
        F[2*n1+0]+=Fe[0]; F[2*n1+1]+=Fe[1]
        F[2*n2+0]+=Fe[2]; F[2*n2+1]+=Fe[3]
    return K, F

def solve_plane_stress_q4_static(env):
    nodes = np.asarray(env['nodes'], float); elems = np.asarray(env['elems'], int)
    E = float(env['E']); nu = float(env['nu'])
    D = (E/(1-nu**2))*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]], float)
    K,F = _assemble_q4_lin(nodes, elems, D, env)
    ndof = K.shape[0]
    fixed_idx, vals = dirichlet_vector(
        env.get('fix', []),
        ndof=ndof,
        ndof_per_node=2,
        dof_map=PLANAR_UV_DOF_MAP,
    )
    free = np.setdiff1d(np.arange(ndof), fixed_idx)
    U = np.zeros(ndof)
    if free.size>0:
        rhs = F[free] - K[np.ix_(free,fixed_idx)] @ vals[fixed_idx]
        U[free] = np.linalg.solve(K[np.ix_(free,free)], rhs)
    U[fixed_idx]=vals[fixed_idx]
    R = K @ U - F
    sigma_node, sigma_elem = _compute_q4_stress(nodes, elems, U, D)
    sx = sigma_node[:,0]; sy = sigma_node[:,1]; txy = sigma_node[:,2]
    sigma_vm = np.sqrt(sx*sx + sy*sy - sx*sy + 3.0*(txy**2))
    return {
        'U':U,
        'R':R,
        'K':K,
        'nodes':nodes,
        'elems':elems,
        'sigma':sigma_node,
        'sigma_elem':sigma_elem,
        'sigma_vm':sigma_vm,
    }

def _compute_q4_stress(nodes, elems, U, D):
    gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
    nn = nodes.shape[0]
    sigma_node = np.zeros((nn,3))
    counts = np.zeros(nn)
    sigma_elem = []
    for conn in elems:
        xe = nodes[conn,:]
        dofs = []
        for a in conn:
            dofs += [2*a, 2*a+1]
        ue = U[dofs]
        sigma_sum = np.zeros(3)
        weight = 0.0
        for xi in gp:
            for eta in gp:
                _, dNdxi = _q4_shape(xi, eta)
                J = dNdxi.T @ xe
                detJ = np.linalg.det(J)
                invJ = np.linalg.inv(J)
                dNdx = dNdxi @ invJ.T
                B = np.zeros((3,8))
                for a in range(4):
                    B[0,2*a] = dNdx[a,0]
                    B[1,2*a+1] = dNdx[a,1]
                    B[2,2*a] = dNdx[a,1]
                    B[2,2*a+1] = dNdx[a,0]
                eps = B @ ue
                sigma_sum += D @ eps * detJ
                weight += detJ
        s_avg = sigma_sum / weight if weight > 0 else np.zeros(3)
        sigma_elem.append(s_avg)
        for a in conn:
            sigma_node[a,:] += s_avg
            counts[a] += 1.0
    counts[counts==0] = 1.0
    sigma_node /= counts[:,None]
    return sigma_node, np.vstack(sigma_elem) if sigma_elem else np.zeros((0,3))
