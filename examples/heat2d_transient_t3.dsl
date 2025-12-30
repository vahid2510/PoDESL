
PROBLEM "Heat : Conduction2D : Transient"
GIVEN
  etype = "T3"
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2],[0,2,3]]
  k = 15.0
  rho = 2700
  c = 900
  t = 0.005
  dt = 0.02
  t_end = 0.2
  q = 1e5
  qedge = [[0,1,2000]]
  fixT = [[0,20]]
  T0 = 20
  store_every = 1
REPORT
  export "heat2d_transient_t3.vtk" vtk
