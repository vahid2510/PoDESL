PROBLEM "Heat : OneD : Linear : Transient"
GIVEN
  nodes = [0.0, 0.5, 1.0]
  elems = [[0,1],[1,2]]
  k = 10.0
  rho = 7800.0
  cp = 500.0
  dt = 0.1
  t_end = 1.0
  T0 = 0.0
  fixT = [[0,100.0]]
REPORT
  print times
