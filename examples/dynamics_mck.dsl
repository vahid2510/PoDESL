PROBLEM "Dynamics : MCKSystem : Transient"
GIVEN
  M = [[2,0],[0,1]]
  C = [[0.02,0],[0,0.02]]
  K = [[2000,-1000],[-1000,1000]]
  u0 = [0,0]
  v0 = [0,0]
  dt = 0.001
  t_end = 0.5
  f = [0.0, 1.0]
REPORT
  print U
