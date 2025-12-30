PROBLEM "Solid : PlaneStress : Q4 : Static"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  bf = [0, -1000]
  fix = [[0,"both",0.0],[3,"both",0.0]]
REPORT
  print U
