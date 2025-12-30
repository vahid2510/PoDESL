PROBLEM "Solid : Truss2D : Static"
GIVEN
  nodes = [[0,0],[1,0],[0.5,0.8]]
  elems = [[0,1],[0,2],[1,2]]
  E = 210e9
  A = 2.0e-4
  supports = [[0,1,1],[1,1,1]]
  loads = [[2, 0.0, -1000.0]]
REPORT
  export "truss2d_linear_U.csv" U
  export "truss2d_linear_forces.csv" element_forces
