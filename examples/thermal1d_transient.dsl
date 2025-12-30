PROBLEM "Thermal : 1DConduction : Transient"
GIVEN
  L = 1.0
  A = 1.0
  k = 12.0
  rho = 7800.0
  c = 500.0
  nel = 8
  dt = 0.1
  t_end = 0.5
  T_left = 100.0
  T_right = 50.0
REPORT
  print T[-1]
