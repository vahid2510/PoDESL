
PROBLEM "Solid : PlaneStress2D : Linear"
GIVEN
  etype = "Q4"
  nodes = [[0,0],[2,0],[2,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.30
  t = 0.01
  traction = [[0, 2, 0.0, 500.0]]
  fix = [[0,'ux'],[0,'uy'],[3,'ux'],[3,'uy']]
REPORT
  export "plate_q4_traction.vtk" vtk
