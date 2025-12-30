
PROBLEM "Solid : PlaneStrain : LinearFE"
GIVEN
  Lx = 48.0e-3
  Ly = 44.0e-3
  t  = 1.0e-3
  nx = 12
  ny = 11
  E  = 70e9
  nu = 0.33
  traction_right_x = 1.0e6
  traction_right_y = 0.0
  clamp_left = True
REPORT
  export "ps_linear_result.vtk" vtk
  export "ps_linear_vm.csv" sigma_vm
