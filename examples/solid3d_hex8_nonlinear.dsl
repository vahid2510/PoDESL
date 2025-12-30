PROBLEM "Solid3D : Hex8 : Nonlinear"
GIVEN
  nodes = [
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
  ]
  elements = [[0,1,2,3,4,5,6,7]]
  material = {"E":210e9,"nu":0.30,"nonlinear_coef":1e6}
  loads = [[6,"uy",-1e5]]
  bcs = [
    {"node":0,"dof":"ux","value":0.0},
    {"node":0,"dof":"uy","value":0.0},
    {"node":0,"dof":"uz","value":0.0},
    {"node":1,"dof":"ux","value":0.0},
    {"node":1,"dof":"uy","value":0.0},
    {"node":1,"dof":"uz","value":0.0},
    {"node":3,"dof":"ux","value":0.0},
    {"node":3,"dof":"uy","value":0.0},
    {"node":3,"dof":"uz","value":0.0},
    {"node":2,"dof":"ux","value":0.0},
    {"node":2,"dof":"uy","value":0.0},
    {"node":2,"dof":"uz","value":0.0}
  ]
  max_iter = 10
  load_steps = 2
REPORT
  print u
