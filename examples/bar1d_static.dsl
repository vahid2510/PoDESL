PROBLEM "Solid : Bar1D : Static"
GIVEN
  L = 2.0
  E = 210e9
  A = 2.0e-4
  nel = 20
  p = 1000.0
  left = "fixed"
  right = "fixed"
REPORT
  print u
  export "bar_u.csv" x u
