PROBLEM "Solid : Plane2D : Q4 : ThermoElastic"
GIVEN
  nodes = [[0,0],[1,0],[1,1],[0,1]]
  elems = [[0,1,2,3]]
  E = 210e9
  nu = 0.3
  alpha = 1.2e-5
  T = [20.0,120.0,120.0,20.0]
  Tref = 20.0
  t = 0.1
  plane = "stress"
  fix = [[0,"both",0.0],[3,"both",0.0]]
REPORT
  print U
