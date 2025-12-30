PROBLEM "Solid : Truss2D : GeomNL"
GIVEN
  nodes = [[0,0],[1,0],[2,0]]
  elems = [[0,1],[1,2]]
  E = 210e9
  A = 1.0e-4
  supports = [[0,1,1],[2,1,1]]
  loads = [[1, 0.0, -2.0e5]]
  CONTROLS = {"method":"arc","load_steps":60,"arc_len":1e-3,"arc_minmax":[1e-5,1e-2],"max_iter":40,"tol":1e-9,"line_search":true}
REPORT
  export "truss2d_snap_path.csv" path
  export "truss2d_snap_U.csv" U
