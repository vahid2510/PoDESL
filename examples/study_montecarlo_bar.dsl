PROBLEM "Study : MonteCarlo : Run"
GIVEN
  child_problem = "Solid : Bar1D : Static"
  base_env = {
    "L": 1.0,
    "A": 1.0e-4,
    "E": 210e9,
    "nel": 4,
    "left": "fixed",
    "right": "traction",
    "traction_right": 1000.0
  }
  param_defs = {
    "traction_right": {"dist": "uniform", "low": 900.0, "high": 1100.0},
    "E": {"dist": "normal", "mean": 210e9, "std": 5e9}
  }
  outputs = ["u"]
  nsamples = 8
  rng_seed = 123
REPORT
  export "study_montecarlo_mean.csv" mean
  export "study_montecarlo_std.csv" std
  print mean
  print std
