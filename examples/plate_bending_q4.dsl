
PROBLEM "Solid : PlateMindlin2D : Linear"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.30
  t = 0.01
  kappa = 0.8333333333
  q = 1e4
  fixW = [[0,0],[1,0],[2,0],[3,0]]
  fixRot = [[0,'rx',0],[0,'ry',0],[1,'rx',0],[1,'ry',0],[2,'rx',0],[2,'ry',0],[3,'rx',0],[3,'ry',0]]
REPORT
  export "plate_bending_q4.vtk" vtk
