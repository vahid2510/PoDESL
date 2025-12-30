PROBLEM "Solid : Truss2D : Static"
GIVEN
  nodes = [[0,0],[3,0],[6,0],[3,3]]
  elems = [
    [0,1, 4.0e-4, 210e9],
    [1,2, 4.0e-4, 210e9],
    [0,3, 4.0e-4, 210e9],
    [1,3, 4.0e-4, 210e9],
    [2,3, 4.0e-4, 210e9]
  ]
  supports = [[0,1,1],[2,1,1]]
  loads = [[3, 0.0, -20000.0]]
REPORT
  print U
  export "truss_U.csv" U
  export "truss_member_forces.csv" member_forces
