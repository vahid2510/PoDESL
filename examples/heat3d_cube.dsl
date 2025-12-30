PROBLEM "Heat : ThreeD : H8 : Transient"
GIVEN
  nodes = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
  elems = [[0,1,2,3,4,5,6,7]]
  k = 10.0
  rho = 7800.0
  cp = 500.0
  dt = 0.1
  t_end = 1.0
  T0 = 0.0
  fixT = [[0,100.0]]
REPORT
  print times
