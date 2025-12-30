
PROBLEM "Materials : VonMises : PlaneStressFE"
GIVEN
  Lx = 1.0
  Ly = 1.0
  t  = 0.01
  nx = 8
  ny = 8
  E = 210e9
  nu = 0.3
  sigma_y0 = 250e6
  H = 1.0e9
  traction_right_x = 0.0
  traction_right_y = -5.0e7
  load_steps = 30
  max_iter   = 40
  tol        = 1e-8
REPORT
  export "ps_result.vtk" vtk
  export "ps_vm.csv" sigma_vm
  export "ps_U.csv" U
