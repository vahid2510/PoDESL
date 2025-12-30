PROBLEM "Materials : VonMises : Uniaxial"
GIVEN
  E = 210e9
  nu = 0.3
  sigma_y0 = 250e6
  H = 1.0e9
  eps_max = 0.02
  nsteps = 200
REPORT
  export "uniaxial_eps.csv" eps
  export "uniaxial_sigma.csv" sigma
  export "uniaxial_epsp.csv" eps_p
  export "uniaxial_Et.csv" Et
