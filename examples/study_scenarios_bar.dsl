PROBLEM "Study : Scenario : Compare"
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
  scenarios = [
    {"name": "baseline", "modifications": []},
    {
      "name": "overload",
      "modifications": [
        {"path": "traction_right", "scale": 1.5}
      ]
    }
  ]
  metrics = [
    {"name": "tip_disp", "expr": "u[-1]"},
    {"name": "max_stress", "expr": "sigma.max()"}
  ]
REPORT
  export "study_scenarios_tip_disp.csv" [s["metrics"]["tip_disp"] for s in scenarios]
  print scenarios
