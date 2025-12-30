
PROBLEM "Solid : Frame2D : PDelta"
GIVEN
  nodes = [[0,0],[3,0],[0,3],[3,3]]
  elems = [[0,1],[0,2],[1,3]]
  E = 210e9
  A = 4.0e-4
  I = 8.0e-6
  udl = [[0, 0.0, 5000.0]]
  point_loads = [[3, 0.0, -30000.0, 0.0]]
  fix = [[0,'ux'],[0,'uy'],[0,'rz'], [1,'ux'],[1,'uy'],[1,'rz']]
  max_iter = 40
  tol = 1.0e-9
REPORT
  print U
