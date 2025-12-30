
PROBLEM "Heat : Conduction2D : Transient"
GIVEN
  etype = "Q4"
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  k = 45.0
  rho = 7800
  c = 500
  t = 0.01
  dt = 0.05
  t_end = 0.5
  theta = 0.5
  fixT = [[0,100],[3,100],[1,0],[2,0]]
  store_every = 1
REPORT
  export "heat2d_transient_q4.vtk" vtk
