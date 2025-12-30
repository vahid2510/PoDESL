PROBLEM "Solid : Plate2D : Mindlin : Modal"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  t = 0.01
  rho = 7800.0
  kappa = 0.8333
  fix = [[0,"all",0.0],[1,"all",0.0],[3,"all",0.0]]
  nmodes = 3
REPORT
  print freq
