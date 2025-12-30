PROBLEM "Solid : Frame2D : Modal : ConsistentMass"
GIVEN
  nodes = [[0,0],[3,0],[3,3]]
  elems = [[0,1],[1,2]]
  E = 210e9
  A = 4e-4
  I = 8e-6
  rho = 7850.0
  fix = [[0,"both",0.0],[1,"uy",0.0]]
  nmodes = 3
REPORT
  print freq
