PROBLEM "Solid : Frame2D : Modal"
GIVEN
  nodes = [[0,0],[3,0]]
  elems = [[0,1]]
  E = 210e9
  A = 4.0e-4
  I = 3.3e-6
  rho = 7850
  supports = [[0,1,1,1]]
REPORT
  export "frame2d_modal_freq.csv" freq
