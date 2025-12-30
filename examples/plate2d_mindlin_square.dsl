PROBLEM "Solid : Plate2D : Mindlin : Static"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  t = 0.01
  p = 1000.0
  fix = [[0,"both",0.0],[1,"both",0.0],[2,"both",0.0],[3,"both",0.0]]
REPORT
  print U
