
PROBLEM "Solid : Elastic3D : Linear"
GIVEN
  etype = "HEX8"
  nodes = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],
           [0,0,1],[1,0,1],[1,1,1],[0,1,1]]
  elems = [[0,1,2,3,4,5,6,7]]
  E = 70e9
  nu = 0.33
  fix = [[0,'ux'],[0,'uy'],[0,'uz'],
         [3,'ux'],[3,'uy'],[3,'uz'],
         [4,'ux'],[4,'uy'],[4,'uz'],
         [7,'ux'],[7,'uy'],[7,'uz']]
  point_loads = [[6, 0, 0, 1e4]]
REPORT
  export "elastic3d_hex8.vtk" vtk
