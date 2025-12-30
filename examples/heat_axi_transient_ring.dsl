PROBLEM "Heat : Axisymmetric : Linear : Transient"
GIVEN
  nodes = [[0.1,0],[0.2,0],[0.2,0.05],[0.1,0.05]]
  elems = [[0,1,2,3]]
  k = 20.0
  rho = 7800.0
  cp = 500.0
  dt = 0.1
  t_end = 1.0
  theta = 1.0
  hedge = [[0,1,10.0,300.0]]
  fix = [[0,300.0],[3,300.0]]
REPORT
  print T_hist
