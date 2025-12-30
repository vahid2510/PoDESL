PROBLEM "Solid3D : Hex8 : Transient"
GIVEN
  nodes = [
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
  ]
  elements = [[0,1,2,3,4,5,6,7]]
  E = 210e9
  nu = 0.30
  rho = 7800.0
  damping_ratio = 0.02
  dt = 1e-4
  t_end = 1e-3
  fixed_dofs = [0,1,2]
  load_func = lambda t: (
    [0.0, 0.0, 0.0, 0.0, 0.0, -5e4*math.sin(2*math.pi*400*t)]
    + [0.0]*18
  )
REPORT
  print result["time"]
  print result["u"][-1]
