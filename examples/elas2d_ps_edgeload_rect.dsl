PROBLEM "Solid : Elasticity2D : PlaneStressEdgeLoad"
GIVEN
  Lx = 2.0
  Ly = 1.0
  nx = 20
  ny = 10
  t = 0.01
  E = 210e9
  nu = 0.3
  clamp_left = 1
  # traction on right edge (downward 5kN/m)
  tx_right = 0.0
  ty_right = -5000.0
REPORT
  export "elas2d_ps_edgeload_U.csv" U
  export "elas2d_ps_edgeload_stress.csv" stress
