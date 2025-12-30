PROBLEM "Solid : Contact1D : Static"
GIVEN
  nn = 3
  springs = [[0,1,1e6],[1,2,1e6]]
  supports = [[0,1]]
  loads = [[2, -5.0e2]]
  contact_gaps = [[1,0, 0.001, 1e9]]
REPORT
  export "contact1d_gap_U.csv" U
  print active_contacts
