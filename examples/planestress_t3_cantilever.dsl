
PROBLEM "Solid : PlaneStress2D : Linear"
GIVEN
  etype = "T3"
  nodes = [[0,0],[2,0],[2,1],[1,1],[0,1]]
  elems = [[0,1,3],[0,3,4],[1,2,3]]
  E = 70e9
  nu = 0.33
  t = 0.02
  b = [0.0, -1000.0]
  fix = [[0,'ux',0],[0,'uy',0],[4,'ux',0],[4,'uy',0]]
REPORT
  print U
