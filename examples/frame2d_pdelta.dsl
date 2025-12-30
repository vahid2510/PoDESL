PROBLEM "Solid : Frame2D : Static"
GIVEN
  nodes = [[0,0],[4,0],[4,3]]
  elems = [[0,1, 3.2e-3, 210e9, 8.0e-6, -2000.0],   # uniform downward local load -2000 N/m
           [1,2, 3.2e-3, 210e9, 8.0e-6,  0.0]]
  supports = [[0,1,1,1],[2,0,1,0]]
  loads = [[1, 0.0, -15000.0, 0.0]]
  nonlinear = true
  max_iter = 50
  tol = 1e-9
REPORT
  print U
  export "frame_member_forces.csv" member_forces
