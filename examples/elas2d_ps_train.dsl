PROBLEM "Solid : Elasticity2D : PlaneStrain"
GIVEN
  Lx = 1.0
  Ly = 1.0
  nx = 8
  ny = 8
  t = 1.0
  E = 70e9
  nu = 0.25
  clamp_left = 1
  loads = [[(nx+1)*ny + nx, 0.0, -1000.0]]
REPORT
  export "elas2d_ps_U.csv" U
  export "elas2d_ps_stress.csv" stress
