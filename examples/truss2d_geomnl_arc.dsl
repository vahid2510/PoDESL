PROBLEM "Solid : Truss2D : GeomNL"
GIVEN
  nodes = [[0,0],[1,0]]
  elems = [[0,1]]
  E = 210e9
  A = 1.0e-4
  supports = [[0,1,1]]
  loads = [[1, 0.0, -5.0e4]]
REPORT
  export "truss2d_geomnl_U.csv" U
  export "truss2d_geomnl_Reac.csv" reactions
