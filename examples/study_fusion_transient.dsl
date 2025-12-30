PROBLEM "Study : Fusion : Transient"
GIVEN
  subsystems = [
    {
      "name": "thermal",
      "solver": "podesl.solvers.heat1d_transient:solve_heat1d_transient",
      "model": {
        "L": 1.0,
        "A": 1.0,
        "k": 10.0,
        "rho": 7800.0,
        "cp": 500.0,
        "nel": 4,
        "nodes": [0.0, 0.25, 0.5, 0.75, 1.0],
        "elems": [[0,1],[1,2],[2,3],[3,4]],
        "dt": 0.1,
        "t_end": 0.0,
        "T_left": 300.0,
        "T_right": 300.0,
        "Tinit": 300.0,
        "fix": [[0, 300.0], [4, 300.0]]
      }
    },
    {
      "name": "bar",
      "solver": "podesl.solvers.bar1d:solve_bar1d_static",
      "model": {
        "L": 1.0,
        "A": 1.0e-4,
        "E": 210e9,
        "nel": 4,
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
      "map": "lambda T: 5e1 * (T[-1][-1] - 300.0)",
      "description": "convert tip temperature into axial load"
    }
  ]
  mode = "transient"
  nsteps = 3
  dt = 0.1
  time_update = "lambda step, t, subsystems: subsystems['thermal'].model.update({'T_right': 300.0 + 5.0 * step})"
REPORT
  export "study_fusion_transient_times.csv" times
  export "study_fusion_transient_tip.csv" [step["u"][-1] for step in history["bar"]]
