PROBLEM "Solid : Plane2D : Q4 : Modal"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  rho = 7800.0
  t = 0.1
  plane = "stress"
  fix = [[0,"both",0.0],[3,"both",0.0]]
  nmodes = 3
REPORT
  print freq
