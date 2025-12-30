PROBLEM "Solid : Axisymmetric : Q4 : Static"
GIVEN
  nodes = [[0.5,0.0],[1.0,0.0],[1.0,1.0],[0.5,1.0]]
  elems = [[0,1,2,3]]
  E = 200e9
  nu = 0.3
  fix = [[0,"ur",0.0],[3,"ur",0.0]]
REPORT
  print U
