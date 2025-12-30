PROBLEM "Solid : Heat2D : Transient"
GIVEN
  Lx=1.0 Ly=1.0 nx=20 ny=20
  kx=45 ky=45 rho=7800 cp=500
  dt=0.02 t_end=0.5 theta=1.0
  T0=100.0
  T_left=0.0
  T_right=0.0
  T_bottom=0.0
  T_top=0.0
REPORT
  export "heat2d_T_final.csv" T
