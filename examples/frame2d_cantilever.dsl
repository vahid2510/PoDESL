PROBLEM "Solid : Frame2D : Static"
GIVEN
  nodes = [[0,0],[2,0]]
  elems = [[0,1]]
  E = 210e9
  A = 4.0e-4
  I = 3.3e-6
  supports = [[0,1,1,1]]
  w_dist = [[0, -1000.0]]
REPORT
  export "frame2d_cantilever_U.csv" U
