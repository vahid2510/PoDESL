PROBLEM "Solid : Elasticity2D : PlaneStress"
GIVEN
  Lx=1.0 Ly=1.0 nx=10 ny=10
  E=70e9 nu=0.33 t=1.0
  alpha=1.2e-5
  dT=50.0
  clamp_left=1
REPORT
  export "elas2d_ps_thermal_U.csv" U
  export "elas2d_ps_thermal_stress.csv" stress
