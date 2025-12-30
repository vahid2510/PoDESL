PROBLEM "Thermal : Heat2D : TransientBC"
GIVEN
  Lx = 1.0
  Ly = 1.0
  nx = 12
  ny = 12
  rho = 7800.0
  c = 500.0
  k = 45.0
  t = 0.01
  q = 0.0
  dt = 0.1
  t_end = 2.0
  theta = 0.5
  T_init = 20.0
  # Dirichlet در چپ
  T_left = 100.0
  # Robin در راست با h(t) خطی بین 10 تا 40
  h_right_series = [[0.0, 10.0], [2.0, 40.0]]
  Tinf_right = 20.0
REPORT
  export "heat2d_transientBC_Tend.csv" T[-1]
