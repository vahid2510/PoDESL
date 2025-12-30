PROBLEM "Frame3D : BeamColumn : Modal"
GIVEN
  nodes = [[0,0,0],[0,0,3]]
  elements = [[0,1]]
  E = 210e9
  G = 80e9
  A = 4e-4
  Iy = 8e-6
  Iz = 8e-6
  J = 1e-5
  rho = 7800.0
  nmodes = 2
  bcs = [
    {"node":0,"dof":"ux","value":0.0},
    {"node":0,"dof":"uy","value":0.0},
    {"node":0,"dof":"uz","value":0.0},
    {"node":0,"dof":"rx","value":0.0},
    {"node":0,"dof":"ry","value":0.0},
    {"node":0,"dof":"rz","value":0.0}
  ]
REPORT
  print freq
