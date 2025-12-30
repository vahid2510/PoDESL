
PROBLEM "Heat : Conduction2D : Linear : Steady"
GIVEN
  etype = "Q4"
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  k = 20.0
  q = 1e5
  hedge = [[0,1,10.0,300.0]]
  fix = [[0,300.0],[3,300.0]]
REPORT
  print T
