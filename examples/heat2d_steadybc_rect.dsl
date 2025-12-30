PROBLEM "Thermal : Heat2D : SteadyBC"
GIVEN
  Lx = 1.0
  Ly = 1.0
  nx = 16
  ny = 16
  k = 15.0
  t = 1.0
  q = 0.0
  # Dirichlet
  T_left = 100.0
  # Convection on right edge
  h_right = 25.0
  Tinf_right = 20.0
  # Insulated top/bottom (هیچ چیزی نده)
REPORT
  export "heat2d_steadybc_T.csv" nodes[:,0] nodes[:,1] T
