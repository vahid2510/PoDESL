PROBLEM "Solid : Truss2D : BucklingLinear"
GIVEN
  nodes = [[0,0],[1,0],[2,0]]
  elems = [[0,1],[1,2]]
  E = 210e9
  A = 1.0e-4
  Nref = -1.0e4
  fix = [[0,"both",0.0],[2,"uy",0.0]]
  nmodes = 3
REPORT
  print lambda
