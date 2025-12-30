PROBLEM "Solid : PlateMindlin : Static"
GIVEN
  Lx = 1.0
  Ly = 1.0
  nx = 12
  ny = 12
  t = 0.02
  E = 210e9
  nu = 0.3
  q = 1000.0     # N/m^2
  clamp_left = 1
  clamp_right = 1
  clamp_bottom = 1
  clamp_top = 1
REPORT
  export "plate_mindlin_U.csv" U
