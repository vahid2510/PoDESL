PROBLEM "Study : Optimize : Gradient"
GIVEN
  child_problem = "Solid : Bar1D : Static"
  base_env = {
    "L": 1.0,
    "A": 1.0e-4,
    "E": 210e9,
    "nel": 4,
    "left": "fixed",
    "right": "traction",
    "traction_right": 500.0
  }
  design_vars = [
    {"path": "traction_right"}
  ]
  x0 = [500.0]
  objective_expr = "(u[-1]*1e6)**2"
  max_iter = 40
  step_size = 0.5
  tol = 1e-8
REPORT
  export "study_optimize_x_opt.csv" x_opt
  export "study_optimize_history.csv" history["f"]
  print x_opt
  print f_opt
