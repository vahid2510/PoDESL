
PROBLEM "Heat : Conduction2D : Steady"
GIVEN
  etype = "Q4"
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  k = 45.0
  t = 0.01
  fixT = [[0,100],[3,100],[1,0],[2,0]]
REPORT
  export "plate_heat_q4.vtk" vtk
