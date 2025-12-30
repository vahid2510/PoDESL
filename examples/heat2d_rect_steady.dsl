PROBLEM "Thermal : Heat2D : Steady"
GIVEN
  Lx = 1.0
  Ly = 1.0
  nx = 10
  ny = 10
  k = 10.0
  t = 1.0
  q = 0.0
  T_left = 100.0
  T_right = 0.0
  T_bottom = 0.0
  T_top = 0.0
REPORT
  export "heat2d_T.csv" nodes[:,0] nodes[:,1] T
