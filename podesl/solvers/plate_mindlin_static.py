import numpy as np
from .mesh_rect import rect_mesh

def solve_plate_mindlin_static(env):
    # DOF per node: [w, thx, thy]
    # Uniform pressure q (positive downward) on area; thickness t, E, nu, shear corr kappa
    if "nodes" in env and "elems" in env:
        nodes = np.array(env["nodes"], dtype=float); elems = np.array(env["elems"], dtype=int)
    else:
        Lx = float(env.get("Lx",1.0)); Ly=float(env.get("Ly",1.0))
        nx = int(env.get("nx",10)); ny=int(env.get("ny",10))
        nodes, elems = rect_mesh(Lx, Ly, nx, ny)
    t = float(env.get("t", 0.01)); E = float(env.get("E", 210e9)); nu = float(env.get("nu", 0.3))
    q = float(env.get("q", 0.0))  # N/m^2 downward
    kappa = float(env.get("kappa", 5.0/6.0))

    # material matrices
    D_b = (E*t**3/(12*(1-nu**2))) * np.array([[1, nu, 0],[nu, 1, 0],[0, 0, (1-nu)/2]])
    G = E/(2*(1+nu)); D_s = kappa * G * t * np.eye(2)

    nn = nodes.shape[0]; ndof = nn*3
    K = np.zeros((ndof, ndof)); F = np.zeros(ndof)

    # 2x2 Gauss for bending; 1-point for shear (to avoid locking)
    a = 1/np.sqrt(3); gp2 = [(-a,-a),(a,-a),(a,a),(-a,a)]
    w2 = [1,1,1,1]

    def shape_Q4(xi,eta):
        N = 0.25*np.array([(1-xi)*(1-eta),
                           (1+xi)*(1-eta),
                           (1+xi)*(1+eta),
                           (1-xi)*(1+eta)], dtype=float)
        dN_dxi = 0.25*np.array([[-(1-eta),  (1-eta), (1+eta), -(1+eta)],
                                [-(1-xi),  -(1+xi), (1+xi),   (1-xi)]], dtype=float)
        return N, dN_dxi

    for e in range(elems.shape[0]):
        conn = elems[e]; xe = nodes[conn,:]
        Ke = np.zeros((12,12)); Fe = np.zeros(12)

        # bending part (2x2)
        for (xi,eta), wt in zip(gp2, w2):
            N, dN_dxi = shape_Q4(xi, eta)
            J = dN_dxi @ xe; detJ = np.linalg.det(J); invJ = np.linalg.inv(J)
            dN_dx = invJ @ dN_dxi  # 2x4
            # curvatures: kx = d(thx)/dx, ky = d(thy)/dy, kxy = d(thx)/dy + d(thy)/dx
            Bb = np.zeros((3,12))
            for i in range(4):
                Bb[0, 3*i+1] = dN_dx[0,i]  # dthx/dx
                Bb[1, 3*i+2] = dN_dx[1,i]  # dthy/dy
                Bb[2, 3*i+1] = dN_dx[1,i]  # dthx/dy
                Bb[2, 3*i+2] = dN_dx[0,i]  # dthy/dx
            Ke += (Bb.T @ D_b @ Bb) * detJ * wt

        # shear part (reduced 1-point at center)
        N, dN_dxi = shape_Q4(0.0, 0.0)
        J = dN_dxi @ xe; detJ = np.linalg.det(J); invJ = np.linalg.inv(J)
        dN_dx = invJ @ dN_dxi
        Bs = np.zeros((2,12))
        for i in range(4):
            # gamma_xz ~ thx + dw/dx ; gamma_yz ~ thy + dw/dy
            Bs[0, 3*i+0] = dN_dx[0,i]   # dw/dx
            Bs[0, 3*i+1] = N[i]         # thx
            Bs[1, 3*i+0] = dN_dx[1,i]   # dw/dy
            Bs[1, 3*i+2] = N[i]         # thy
        Ke += (Bs.T @ D_s @ Bs) * detJ

        # uniform pressure load q => Fe_w = \int N q dΩ
        if abs(q) > 0:
            for (xi,eta), wt in zip(gp2, w2):
                N, dN_dxi = shape_Q4(xi, eta)
                J = dN_dxi @ xe; detJ = np.linalg.det(J)
                for i in range(4):
                    Fe[3*i+0] += N[i] * q * detJ * wt

        # assemble
        dofs = []
        for n in conn: dofs.extend([3*n+0,3*n+1,3*n+2])
        for i in range(12):
            for j in range(12):
                K[dofs[i], dofs[j]] += Ke[i,j]
        for i in range(12): F[dofs[i]] += Fe[i]

    # supports
    fixed = []
    if "supports" in env:
        for s in np.array(env["supports"], dtype=float):
            n = int(s[0]); fix_w=int(s[1]); fix_tx=int(s[2]); fix_ty=int(s[3])
            if fix_w: fixed.append(3*n+0)
            if fix_tx: fixed.append(3*n+1)
            if fix_ty: fixed.append(3*n+2)

    # convenience clamps on edges
    xs=nodes[:,0]; ys=nodes[:,1]; tol=1e-12
    minx,maxx,miny,maxy = xs.min(),xs.max(),ys.min(),ys.max()
    edges = {
        "clamp_left":  np.where(np.abs(xs-minx)<tol)[0],
        "clamp_right": np.where(np.abs(xs-maxx)<tol)[0],
        "clamp_bottom":np.where(np.abs(ys-miny)<tol)[0],
        "clamp_top":   np.where(np.abs(ys-maxy)<tol)[0],
    }
    for name,ids in edges.items():
        if int(env.get(name,0))==1:
            for i in ids:
                fixed.extend([3*i+0,3*i+1,3*i+2])

    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed)
    if free.size==0: raise ValueError("No free DOFs. Check clamps/supports.")
    U = np.zeros(ndof); U[free] = np.linalg.solve(K[free][:,free], F[free])
    return {"nodes": nodes, "elems": elems, "U": U}
