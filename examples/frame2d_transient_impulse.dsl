PROBLEM "Solid : Frame2D : Transient"
GIVEN
  nodes = [[0,0],[3,0]]
  elems = [[0,1]]
  E = 210e9
  A = 4.0e-4
  I = 3.3e-6
  rho = 7850
  supports = [[0,1,1,1],[1,0,1,0]]
  dt=0.002 t_end=0.2 beta=0.25 gamma=0.5
  rayleigh_a0=0.0 rayleigh_a1=0.0005
  loads_t = [[0.0, 1, 0, -1e3, 0],
             [0.02,1, 0, -1e3, 0],
             [0.021,1, 0, 0,    0],
             [0.2, 1, 0, 0,    0]]
REPORT
  export "frame2d_transient_Uhist.npy" U_hist
