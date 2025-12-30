PROBLEM "Solid : Truss3D : Static"
GIVEN
  nodes = [[0,0,0],[1,0,0],[1,1,0]]
  elems = [[0,1],[1,2]]
  E = 210e9
  A = 1.0e-4
  loads = [[2, 0.0, -1000.0, 0.0]]
  fix = [[0,"both",0.0],[1,"uy",0.0]]
REPORT
  print U
