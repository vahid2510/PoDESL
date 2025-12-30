PROBLEM "Solid : Elasticity2D : PlaneStrain"
GIVEN
  Lx=1.0 Ly=0.5 nx=12 ny=6
  E=210e9 nu=0.3 t=1.0
  bx=0.0 by=-1000.0
  clamp_left=1
  tx_right=1e5
  ty_right=0.0
REPORT
  export "elas2d_plane_strain_U.csv" U
  export "elas2d_plane_strain_stress.csv" stress
