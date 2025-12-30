
PROBLEM "Solid : Truss2D : Static"
GIVEN
  nodes = [[0,0],[1,0],[2,0],[3,0],[4,0]]
  elems = [[0,1],[1,2],[2,3],[3,4]]
  E = 210e9
  A = 3.0e-4
  point_loads = [[4, 0.0, -1000.0]]
  fix = [[0,'x'],[0,'y']]
REPORT
  export "truss2d_cantilever.vtk" vtk
