
PROBLEM "Solid : PlaneStress2D : Linear"
GIVEN
  etype = "Q4"
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  t = 0.01
  tedge = [[0,1,0.0,-1e6]]
  fix = [[0,'ux',0.0],[0,'uy',0.0],[3,'ux',0.0],[3,'uy',0.0]]
REPORT
  print U
