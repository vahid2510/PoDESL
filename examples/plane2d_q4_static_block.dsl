PROBLEM "Solid : Plane2D : Q4 : Static"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  t = 0.1
  plane = "stress"
  fix = [[0,"both",0.0],[3,"both",0.0]]
  loads = [[1,1000.0,0.0],[2,1000.0,0.0]]
REPORT
  print U
