
PROBLEM "Solid : PlaneStrain2D : Linear"
GIVEN
  etype = "T3"
  nodes = [[0,0],[2,0],[2,1],[0,1]]
  elems = [[0,1,2],[0,2,3]]
  E = 70e9
  nu = 0.33
  t = 0.05
  traction = [[0, 1, 0.0, 5e5]]
  fix = [[0,'ux'],[0,'uy'],[3,'ux'],[3,'uy']]
REPORT
  export "plane_strain_t3.vtk" vtk
