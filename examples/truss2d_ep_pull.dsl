PROBLEM "Solid : Truss2D : ElastoPlastic"
GIVEN
  nodes = [[0,0],[1,0]]
  elems = [[0,1]]
  E = 210e9
  A = 1.0e-4
  sigy = 250e6
  H = 1.0e9
  supports = [[0,1,1]]
  loads = [[1, 1.0e5, 0.0]]
REPORT
  export "truss2d_ep_U.csv" U
  export "truss2d_ep_N.csv" N
  export "truss2d_ep_ep.csv" ep
