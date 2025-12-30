PROBLEM "Heat : Planar : Linear : Transient"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  k = 20.0
  rho = 7800.0
  cp = 500.0
  dt = 0.1
  t_end = 1.0
  theta = 1.0
  hedge = [[0,1,10.0,300.0]]
  fix = [[3,400.0]]
REPORT
  print T_hist
