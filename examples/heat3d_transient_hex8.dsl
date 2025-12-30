
PROBLEM "Heat : Conduction3D : Transient"
GIVEN
  etype = "HEX8"
  nodes = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
  elems = [[0,1,2,3,4,5,6,7]]
  k = 45.0
  rho = 7800
  c = 500
  dt = 0.05
  t_end = 0.5
  fixT = [[0,100],[3,100],[4,100],[7,100],[1,0],[2,0],[5,0],[6,0]]
  store_every = 1
REPORT
  export "heat3d_transient_hex8.vtk" vtk
