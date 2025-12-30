import numpy as np

from podesl.uncertainty import mc_propagate
from podesl.optimize import optimize_design, finite_diff_gradient
from podesl.scenario import apply_modifications, run_scenarios
from podesl.safety import check_assertions, compute_safety_score
from podesl.calibrate import calibrate_parameters
from podesl.units import Quantity, strip_units
from podesl.adapt import adapt_loop


def test_mc_propagate_stats():
    def solver(model):
        a = model["a"]
        b = model["b"]
        return {"sum": a + b}

    base = {"a": 1.0, "b": 2.0}
    params = {"a": {"dist": "normal", "mean": 1.0, "std": 0.0}}
    stats = mc_propagate(solver, base, params, outputs=["sum"], nsamples=4, rng_seed=42)
    assert stats["samples"] == 4
    assert np.allclose(stats["mean"]["sum"], 3.0)


def test_optimize_design_quadratic():
    f = lambda x: float((x[0] - 1.0) ** 2)
    res = optimize_design(f, np.array([5.0]), step_size=0.1, max_iter=50, tol=1e-8)
    assert abs(res["x_opt"][0] - 1.0) < 1e-3


def test_finite_diff_gradient():
    f = lambda x: float(x[0] ** 2 + 2.0 * x[1])
    grad = finite_diff_gradient(f, np.array([3.0, -4.0]))
    assert np.allclose(grad, np.array([6.0, 2.0]), atol=1e-5)


def test_scenario_manager():
    base = {"loads": [{"value": 10.0}]}
    scenarios = [
        {"name": "scale", "modifications": [{"path": "loads[0].value", "scale": 2.0}]}
    ]

    def solver(model):
        return {"loads": model["loads"], "sum": model["loads"][0]["value"]}

    metrics = {"sum": lambda res: res["sum"]}
    summary = run_scenarios(base, scenarios, solver, metrics)
    assert summary["scenarios"][0]["metrics"]["sum"] == 20.0

def test_apply_modifications_set_add_scale():
    base = {"value": 1.0, "loads": [{"value": 10.0}]}
    mods = [
        {"path": "value", "set": 2.0},
        {"path": "loads[0].value", "add": 5.0},
        {"path": "loads[0].value", "scale": 2.0},
    ]
    updated = apply_modifications(base, mods)
    assert updated["value"] == 2.0
    assert updated["loads"][0]["value"] == ( (10.0 + 5.0) * 2.0 )


def test_safety_assertions():
    result = {"stress": np.array([10.0, 5.0]), "freq": np.array([22.0])}
    assertions = [
        {"type": "max_less_than", "field": "stress", "limit": 12.0},
        {"type": "value_in_range", "expr": "freq[0]", "min": 20.0, "max": 25.0},
    ]
    reports = check_assertions(result, assertions)
    assert all(r["ok"] for r in reports)
    assert compute_safety_score(reports) == 1.0


def test_calibrate_parameters_linear():
    data_x = np.linspace(0.0, 1.0, 5)
    data_y = 3.0 * data_x

    def model_builder(params):
        return {"k": params["k"]}

    def solver(model):
        return model

    def obs_extractor(res, x):
        return res["k"] * x

    res = calibrate_parameters(
        ["k"],
        np.array([1.0]),
        model_builder,
        solver,
        data_x,
        data_y,
        obs_extractor,
        max_iter=200,
        step_size=0.1,
    )
    assert abs(res["param_opt"][0] - 3.0) < 1e-2


def test_units_and_strip():
    q1 = Quantity(2.0, "m")
    q2 = Quantity(3.0, "m")
    q3 = q1 + q2
    assert q3.value == 5.0
    q4 = q1 * Quantity(4.0, "N")
    assert "m" in q4.unit and "N" in q4.unit
    data = {"value": q3, "list": [q1, 5.0]}
    stripped = strip_units(data)
    assert stripped["value"] == 5.0 and stripped["list"][0] == 2.0


def test_adapt_loop_stops():
    def solver(model):
        return {"u": model["mesh_size"]}

    def functional(res):
        return res["u"]

    def estimator(model, result, func):
        return model["mesh_size"]

    def refine(model, result, err):
        new = model.copy()
        new["mesh_size"] *= 0.5
        return new

    history = adapt_loop({"mesh_size": 0.2}, solver, estimator, refine, functional, max_refine=5, tol=0.05)
    assert history["final_model"]["mesh_size"] <= 0.05
