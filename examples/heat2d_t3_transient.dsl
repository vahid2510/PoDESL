
PROBLEM "Heat : Conduction2D : Linear : Transient"
GIVEN
  etype = "T3"
  nodes = [[0,0],[2,0],[2,1],[0,1]]
  elems = [[0,1,2],[0,2,3]]
  k = 15.0
  rho = 8900.0
  cp = 385.0
  dt = 0.05
  t_end = 0.5
  theta = 1.0
  qedge = [[0,1,1000.0]]
  fix = [[0,300.0],[3,300.0]]
REPORT
  print T_hist
