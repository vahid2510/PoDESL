PROBLEM "Solid : PlaneStrain : Q4 : Static"
GIVEN
  nodes = [[0,0],[2,0],[2,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 30e9
  nu = 0.25
  fix = [[0,"both",0.0],[3,"both",0.0]]
REPORT
  print U
