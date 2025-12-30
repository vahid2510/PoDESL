
PROBLEM "Solid : Truss3D : Linear"
GIVEN
  nodes = [[0,0,0],[1,0,0]]
  elems = [[0,1]]
  A = 1.0e-4
  E = 200e9
  fix = [[0,'ux'],[0,'uy'],[0,'uz']]
  point_loads = [[1, 1000.0, 0.0, 0.0]]
REPORT
  print U
