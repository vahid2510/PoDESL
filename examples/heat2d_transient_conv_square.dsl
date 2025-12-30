PROBLEM "Solid : Heat2D : TransientConv"
GIVEN
  Lx=1.0 Ly=1.0 nx=20 ny=20
  kx=45 ky=45 rho=7800 cp=500
  dt=0.02 t_end=0.5 theta=1.0
  T0=100.0
  qdot=2.0e5
  h_left=50  Tinf_left=20
  h_right=50 Tinf_right=20
  h_bottom=50 Tinf_bottom=20
  h_top=50   Tinf_top=20
REPORT
  export "heat2d_conv_Thist.npy" T_hist
