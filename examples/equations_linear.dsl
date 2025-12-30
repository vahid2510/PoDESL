PROBLEM "LinearAlgebra : Equations : Direct"
EQUATIONS
  unknowns = ["u0","u1","u2"]
  eq " 2*u0 - 1*u1      =  5"
  eq "-1*u0 + 2*u1 - 1*u2 = -1"
  eq "        -1*u1 + 2*u2 =  1"
REPORT
  print x
  export "eqs_x.csv" x
