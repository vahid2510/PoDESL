
PROBLEM "Heat : Conduction2D : Steady"
GIVEN
  etype = "T3"
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2],[0,2,3]]
  k = 20.0
  t = 0.02
  qedge = [[0, 1, 1000.0]]
  fixT = [[0,0],[1,0]]
REPORT
  export "plate_heat_t3.vtk" vtk
