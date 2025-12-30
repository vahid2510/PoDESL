PROBLEM "Solid : Frame2D : StaticNL"
GIVEN
  nodes = [[0,0],[3,0],[3,3]]
  elems = [[0,1],[1,2]]
  E = 210e9
  A = 4.0e-4
  I = 3.3e-6
  supports = [[0,1,1,1],[2,0,1,0]]
  loads = [[1, 0.0, -2.0e4, 0.0]]
  CONTROLS = {"method":"loadstep","load_steps":20,"max_iter":30,"tol":1e-9,"line_search":true}
REPORT
  export "frame2d_loadstep_path.csv" path
  export "frame2d_loadstep_U.csv" U
