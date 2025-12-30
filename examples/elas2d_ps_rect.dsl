PROBLEM "Solid : Elasticity2D : PlaneStress"
GIVEN
  Lx = 2.0
  Ly = 1.0
  nx = 12
  ny = 6
  t = 0.01
  E = 210e9
  nu = 0.3
  clamp_left = 1
  supports = []
  loads = [[(nx+1)*ny + nx, 0.0, -1000.0]]
REPORT
  export "elas2d_U.csv" U
  export "elas2d_stress.csv" stress
