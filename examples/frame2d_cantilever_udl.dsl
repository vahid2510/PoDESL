PROBLEM "Solid : Frame2D : Static"
GIVEN
  nodes = [[0,0],[4,0]]
  elems = [[0,1]]
  E = 210e9
  A = 4e-4
  I = 8e-6
  edist = [[0, 0.0, -1000.0]]
  fix = [[0,"both",0.0]]
REPORT
  print U
