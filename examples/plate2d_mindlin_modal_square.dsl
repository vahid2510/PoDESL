PROBLEM "Solid : Plate2D : Mindlin : Modal"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  t = 0.01
  rho = 7800.0
  kappa = 0.8333333
  fix = [[0,"both",0.0],[1,"both",0.0],[2,"both",0.0],[3,"both",0.0]]
  nmodes = 6
REPORT
  print freq
