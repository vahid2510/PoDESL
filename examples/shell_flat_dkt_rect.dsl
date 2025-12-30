PROBLEM "Solid : Shell : FlatDKT"
GIVEN
  Lx = 1.0
  Ly = 1.0
  nx = 10
  ny = 10
  t = 0.01
  E = 70e9
  nu = 0.25
  q = 500.0
  clamp_left = 1
  clamp_right = 1
  clamp_bottom = 1
  clamp_top = 1
REPORT
  export "shell_flat_dkt_U.csv" U
