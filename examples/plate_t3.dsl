
PROBLEM "Solid : PlaneStress2D : Linear"
GIVEN
  etype = "T3"
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2],[0,2,3]]
  E = 210e9
  nu = 0.30
  t = 0.01
  point_loads = [[2, 0.0, -1000.0]]
  fix = [[0,'ux'],[0,'uy'],[3,'ux'],[3,'uy']]
REPORT
  export "plate_t3.vtk" vtk
