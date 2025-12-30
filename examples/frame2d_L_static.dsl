PROBLEM "Solid : Frame2D : Static"
GIVEN
  nodes = [[0,0],[3,0],[3,3]]
  elems = [[0,1],[1,2]]
  E = 210e9
  A = 4e-4
  I = 8e-6
  loads = [[2, 0.0, -10000.0, 0.0]]
  fix = [[0, "both", 0.0], [1, "uy", 0.0]]
REPORT
  print U
