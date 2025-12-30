PROBLEM "Solid : Frame2D : Buckling : Linear"
GIVEN
  nodes = [[0,0],[3,0],[3,3]]
  elems = [[0,1],[1,2]]
  E = 210e9
  A = 4e-4
  I = 8e-6
  Nref = [1.0e5, 2.0e5]
  fix = [[0,"both",0.0],[1,"uy",0.0]]
  nmodes = 3
REPORT
  print lambda_cr
