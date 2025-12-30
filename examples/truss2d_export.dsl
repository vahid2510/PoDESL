PROBLEM "Solid : Truss2D : Static"
GIVEN
  nodes = [[0,0],[1,0],[1,1]]
  elems = [[0,1, 1e-4, 210e9],
           [1,2, 1e-4, 210e9],
           [0,2, 1e-4, 210e9]]
  supports = [[0, 1,1], [1, 0,1]]
  loads = [[2, 0, -1000]]
REPORT
  export "truss_u.csv" U
