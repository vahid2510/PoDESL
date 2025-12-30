PROBLEM "Solid : Frame2D : StaticNL"
GIVEN
  nodes = [[0,0],[0,3]]
  elems = [[0,1]]
  E = 210e9
  A = 6.0e-4
  I = 2.0e-6
  supports = [[0,1,1,1],[1,0,1,0]]
  loads = [[1, 0.0, -1.0e4, 0.0]]
  CONTROLS = {"method":"arc","load_steps":50,"arc_len":5e-4,"max_iter":30,"tol":1e-8,"line_search":true}
REPORT
  export "frame2d_postbuckling_path.csv" path
  export "frame2d_postbuckling_U.csv" U
