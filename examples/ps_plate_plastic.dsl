PROBLEM "Materials : VonMises : PlaneStressFE"
GIVEN
  Lx = 1.0
  Ly = 1.0
  t = 0.01
  nx = 1
  ny = 1
  E = 210e9
  nu = 0.3
  sigma_y0 = 250e6
  H = 1.0e9
  ux_bar = 0.004
  load_steps = 40
  max_iter = 40
  tol = 1e-8
  penalty = 1e9
REPORT
  export "ps_path.csv" path
  export "ps_U.csv" U
