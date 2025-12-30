PROBLEM "Solid : Beam : Modal"
GIVEN
  L = 2.0
  E = 210e9
  I = 8.0e-6
  rho = 7800.0
  A = 3.0e-4
  nel = 40
  left = "clamped"
  right = "free"
  nmodes = 4
REPORT
  print freq
