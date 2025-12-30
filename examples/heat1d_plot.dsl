PROBLEM "Thermal : Heat1D : Steady"
GIVEN
  L = 1.0
  k = 200.0
  A = 1.0
  qdot = 0.0
  nel = 20
  left = "Dirichlet"; T_left = 100.0
  right = "Robin"; h_right = 10.0; Tinf_right = 20.0
REPORT
  plot x T
  export "heat_T.csv" x T
