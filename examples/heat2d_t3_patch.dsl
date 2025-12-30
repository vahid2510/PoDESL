
PROBLEM "Heat : Conduction2D : Linear : Steady"
GIVEN
  etype = "T3"
  nodes = [[0,0],[2,0],[2,1],[0,1]]
  elems = [[0,1,2],[0,2,3]]
  k = 45.0
  qedge = [[0,1,1000.0]]
  fix = [[0,350.0],[3,350.0]]
REPORT
  print T
