
PROBLEM "Solid : Elastic3D : Linear"
GIVEN
  etype = "TET4"
  nodes = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
  elems = [[0,1,2,3]]
  E = 200e9
  nu = 0.28
  body = [0,0,-9.81]
  rho = 7800
  fix = [[0,'ux'],[0,'uy'],[0,'uz']]
REPORT
  export "elastic3d_tet4.vtk" vtk
