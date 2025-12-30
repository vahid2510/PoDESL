
PROBLEM "Solid : PlaneStrain : LinearFE"
GIVEN
  Lx = 48.0e-3
  Ly = 44.0e-3
  t  = 1.0e-3
  nx = 24
  ny = 22
  E  = 70e9
  nu = 0.3
  traction_right_x = 1.0e6
  traction_right_y = 0.0
  clamp_left = True
REPORT
  export "cook_ps_linear.vtk" vtk
