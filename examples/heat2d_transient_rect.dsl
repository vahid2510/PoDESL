PROBLEM "Thermal : Heat2D : Transient"
GIVEN
  Lx = 1.0
  Ly = 1.0
  nx = 10
  ny = 10
  rho = 7800.0
  c = 500.0
  k = 45.0
  t = 0.01
  q = 0.0
  dt = 0.05
  t_end = 1.0
  theta = 0.5
  T_init = 20.0
  T_left = 100.0
  T_right = 20.0
REPORT
  export "heat2d_transient_Tend.csv" T[-1]
