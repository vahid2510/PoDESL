PROBLEM "Solid : Truss2D : Static"
GIVEN
  nodes = [[0,0],[1,0],[0,1]]
  elems = [[0,1],[1,2],[2,0]]
  E = [210e9, 210e9, 210e9]
  A = [8e-5, 8e-5, 8e-5]
  loads = [[1, 1000.0, -500.0]]
  fix = [[0,"both",0.0],[2,"uy",0.0]]
REPORT
  print U
  print element_forces
