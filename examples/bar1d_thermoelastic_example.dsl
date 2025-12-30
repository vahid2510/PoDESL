PROBLEM "Solid : Bar1D : ThermoElastic : Static"
GIVEN
  nodes = [0.0, 1.0]
  elems = [[0,1]]
  E = 210e9
  A = 1.0e-4
  alpha = 1.2e-5
  T = [20.0, 120.0]
  Tref = 20.0
  fix = [[0,0.0]]
REPORT
  print U
