PROBLEM "Study : Fusion : Steady"
GIVEN
  subsystems = [
    {
      "name": "thermal",
      "solver": "podesl.solvers.heat1d:solve_heat1d_steady",
      "model": {
        "L": 1.0,
        "A": 1.0,
        "k": 10.0,
        "nel": 10,
        "T_left": 300.0,
        "T_right": 350.0
      }
    },
    {
      "name": "bar",
      "solver": "podesl.solvers.bar1d:solve_bar1d_static",
      "model": {
        "L": 1.0,
        "A": 1.0e-4,
        "E": 210e9,
        "nel": 10,
        "left": "fixed",
        "right": "traction",
        "traction_right": 0.0
      }
    }
  ]
  couplings = [
    {
      "src": "thermal",
      "src_field": "T",
      "dst": "bar",
      "dst_field": "traction_right",
      "map": "lambda T: 1.0e2 * (T[-1] - 300.0)",
      "description": "convert temp difference to axial load"
    }
  ]
  mode = "steady"
  max_iter = 3
  tol = 1e-6
REPORT
  print results["bar"]["u"][-1]
