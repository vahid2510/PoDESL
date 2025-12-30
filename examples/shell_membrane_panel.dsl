PROBLEM "Shell : Membrane : Static"
GIVEN
  nodes = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0]
  ]
  elements = [
    [0, 1, 2]
  ]
  thickness = 0.01
  E = 210e9
  nu = 0.3
  bcs = [
    {"node": 0, "dof": "ux", "value": 0.0},
    {"node": 0, "dof": "uy", "value": 0.0},
    {"node": 1, "dof": "uy", "value": 0.0}
  ]
  loads = [
    {"node": 2, "dof": "uy", "value": -1000.0}
  ]
REPORT
  print u
  print element_stress
