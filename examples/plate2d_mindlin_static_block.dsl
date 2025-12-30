PROBLEM "Solid : Plate2D : Mindlin : Static"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  t = 0.01
  kappa = 0.8333
  fix = [[0,"all",0.0],[1,"all",0.0],[3,"all",0.0]]
  loads = [[2,-1000.0]]
REPORT
  print U
