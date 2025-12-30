
PROBLEM "Heat : Conduction3D : Transient"
GIVEN
  etype = "TET4"
  nodes = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
  elems = [[0,1,2,3]]
  k = 15.0
  rho = 2700
  c = 900
  dt = 0.02
  t_end = 0.2
  q = 1e5
  qface = [[0,0,2000]]
  fixT = [[0,20]]
  T0 = 20
  store_every = 1
REPORT
  export "heat3d_transient_tet4.vtk" vtk
