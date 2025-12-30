import numpy as np

from podesl.solvers.study_drivers import (
    solve_study_montecarlo,
    solve_study_optimize,
    solve_study_scenarios,
)

BASE_BAR = {
    "L": 1.0,
    "A": 1.0e-4,
    "E": 210e9,
    "nel": 4,
    "left": "fixed",
    "right": "traction",
    "traction_right": 1000.0,
}


def _bar_env():
    return {key: value for key, value in BASE_BAR.items()}


def test_study_montecarlo_wrapper():
    env = {
        "child_problem": "Solid : Bar1D : Static",
        "base_env": _bar_env(),
        "param_defs": {
            "traction_right": {"dist": "normal", "mean": 1000.0, "std": 0.0},
            "E": {"dist": "normal", "mean": 210e9, "std": 0.0},
        },
        "outputs": ["u"],
        "nsamples": 3,
    }
    stats = solve_study_montecarlo(env)
    assert stats["samples"] == 3
    tip_mean = stats["mean"]["u"][-1]
    assert np.isclose(tip_mean, stats["raw"]["u"][0, -1])


def test_study_optimize_wrapper():
    env = {
        "child_problem": "Solid : Bar1D : Static",
        "base_env": _bar_env(),
        "design_vars": ["traction_right"],
        "x0": [800.0],
        "objective_expr": "(u[-1]*1e6)**2",
        "step_size": 0.5,
        "max_iter": 40,
        "tol": 1e-8,
    }
    res = solve_study_optimize(env)
    initial_tip = env["x0"][0] / (BASE_BAR["E"] * BASE_BAR["A"])
    optimized_tip = res["result"]["u"][-1]
    assert abs(optimized_tip) < abs(initial_tip)
    assert abs(res["x_opt"][0]) < env["x0"][0]


def test_study_scenarios_wrapper():
    env = {
        "child_problem": "Solid : Bar1D : Static",
        "base_env": _bar_env(),
        "scenarios": [
            {"name": "baseline", "modifications": []},
            {
                "name": "boost",
                "modifications": [{"path": "traction_right", "scale": 2.0}],
            },
        ],
        "metrics": [
            {"name": "tip", "expr": "u[-1]"},
            {"name": "max_sigma", "expr": "sigma.max()"},
        ],
    }
    summary = solve_study_scenarios(env)
    assert len(summary["scenarios"]) == 2
    tip_base = summary["scenarios"][0]["metrics"]["tip"]
    tip_boost = summary["scenarios"][1]["metrics"]["tip"]
    assert tip_boost > tip_base
