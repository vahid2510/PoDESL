PROBLEM "Solid : Shell2D : Static"
GIVEN
  Lx=1.0 Ly=1.0 nx=8 ny=8
  E=210e9 nu=0.3 t=0.01
  q=1000.0
  clamp_left=1
  clamp_right=1
  clamp_bottom=1
  clamp_top=1
REPORT
  export "shell2d_plate_U.csv" U
