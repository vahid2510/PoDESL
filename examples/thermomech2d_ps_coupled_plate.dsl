PROBLEM "Solid : ThermoMech2D : PS"
GIVEN
  Lx=1.0 Ly=1.0 nx=12 ny=12
  E=70e9 nu=0.33 t=1.0
  alpha=1.2e-5
  dT=60.0
  supports = [[0,1,1]]
REPORT
  export "thermomech2d_ps_U.csv" U
