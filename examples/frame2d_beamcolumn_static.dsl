PROBLEM "Frame2D : BeamColumn : Static"
GIVEN
  nodes = [[0,0],[0,3]]
  elements = [[0,1]]
  E = 210e9
  A = 4e-4
  I = 8e-6
  fix = [[0,"ux",0.0],[0,"uy",0.0],[0,"rz",0.0]]
  loads = [[1,0.0,-1000.0,0.0]]
REPORT
  print U
