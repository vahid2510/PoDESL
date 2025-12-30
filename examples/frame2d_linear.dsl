PROBLEM "Solid : Frame2D : Static"
GIVEN
  nodes = [[0,0],[4,0],[4,3]]
  # elems: n1, n2, A, E, I, wloc(optional, local -y). Here no distributed load.
  elems = [[0,1, 3.2e-3, 210e9, 8.0e-6],
           [1,2, 3.2e-3, 210e9, 8.0e-6]]
  # supports: node, fix_u, fix_v, fix_theta
  supports = [[0,1,1,1],[2,0,1,0]]
  # nodal loads: node, Fx, Fy, Mz
  loads = [[1, 0.0, -10000.0, 0.0]]
REPORT
  print U
  export "frame_U.csv" U
