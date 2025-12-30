import numpy as np
def rect_mesh(Lx, Ly, nx, ny):
    nx=int(nx); ny=int(ny)
    xs = np.linspace(0.0, Lx, nx+1)
    ys = np.linspace(0.0, Ly, ny+1)
    nodes = np.array([[x,y] for j in range(ny+1) for x in xs for y in [ys[j]]], dtype=float)
    def nid(i,j): return j*(nx+1)+i
    elems = []
    for j in range(ny):
        for i in range(nx):
            n00=nid(i,j); n10=nid(i+1,j); n11=nid(i+1,j+1); n01=nid(i,j+1)
            elems.append([n00,n10,n11,n01])
    return nodes, np.array(elems,dtype=int)
