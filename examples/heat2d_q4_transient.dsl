
PROBLEM "Heat : Conduction2D : Linear : Transient"
GIVEN
  etype = "Q4"
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  k = 35.0
  rho = 7800.0
  cp = 500.0
  dt = 0.1
  t_end = 1.0
  theta = 0.5
  hedge = [[0,1,20.0,300.0]]
  fix = [[0,300.0],[3,300.0]]
REPORT
  print T_hist
