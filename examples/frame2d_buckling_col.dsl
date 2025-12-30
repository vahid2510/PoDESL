PROBLEM "Solid : Frame2D : Buckling"
GIVEN
  nodes = [[0,0],[3,0]]
  elems = [[0,1]]
  E = 210e9
  A = 4.0e-4
  I = 3.3e-6
  supports = [[0,1,1,1],[1,0,1,0]]
  N_axial = [1.0e5]
REPORT
  export "frame2d_buckling_lambda.csv" lambda_cr
