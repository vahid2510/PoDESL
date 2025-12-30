PROBLEM "Solid : Truss2D : Nonlinear"
GIVEN
  nodes = [[0,0],[1,0]]
  elems = [[0,1]]
  E = 210e9
  A = 2.0e-4
  supports = [[0,1,1]]
  loads = [[1, 10000.0, 0.0]]
  max_iter = 40
  tol = 1e-10
REPORT
  export "truss2d_nonlinear_U.csv" U
  export "truss2d_nonlinear_N.csv" N
