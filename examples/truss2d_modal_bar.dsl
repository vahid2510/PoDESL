PROBLEM "Solid : Truss2D : Modal"
GIVEN
  nodes = [[0,0],[2,0]]
  elems = [[0,1]]
  E=210e9 A=1.0e-4 rho=7800
  supports = [[0,1,1]]
REPORT
  export "truss2d_modal_freq.csv" freq
