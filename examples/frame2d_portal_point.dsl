
PROBLEM "Solid : Frame2D : Linear"
GIVEN
  nodes = [[0,0],[0,3],[4,3]]
  elems = [[0,1],[1,2]]
  E = 200e9
  A = 1.0e-3
  I = 2.0e-6
  point_loads = [[2, 0.0, -1000.0, 0.0]]
  fix = [[0,'ux'],[0,'uy'],[0,'rz']]
REPORT
  print U
