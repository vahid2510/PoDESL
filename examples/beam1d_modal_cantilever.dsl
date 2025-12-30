PROBLEM "Solid : Beam1D : Modal"
GIVEN
  nodes = [0.0, 0.5, 1.0]
  elems = [[0,1],[1,2]]
  E = 210e9
  I = 8.0e-6
  rho = 7800.0
  A = 4.0e-4
  fix = [[0,"both",0.0]]
  nmodes = 3
REPORT
  print freq
