import numpy as np

from podesl.solvers.study_fusion import solve_study_fusion


def _thermal_solver(model):
    temp = model.get("temperature", 300.0)
    return {"temperature": temp}


def _solid_solver(model):
    load = model.get("load", 0.0)
    return {"strain": load * 0.1}


def test_study_fusion_steady_custom_callables():
    env = {
        "subsystems": [
            {"name": "thermal", "solver": _thermal_solver, "model": {"temperature": 320.0}},
            {"name": "solid", "solver": _solid_solver, "model": {"load": 0.0}},
        ],
        "couplings": [
            {
                "src": "thermal",
                "src_field": "temperature",
                "dst": "solid",
                "dst_field": "load",
                "map": "lambda temp: scale * (temp - 300.0)",
            }
        ],
        "context": {"scale": 2.0},
        "mode": "steady",
        "max_iter": 2,
    }
    result = solve_study_fusion(env)
    assert result["mode"] == "steady"
    solid_strain = result["results"]["solid"]["strain"]
    assert np.isclose(solid_strain, 40.0 * 0.1)
    assert len(result["residual_history"]) >= 1


def test_study_fusion_transient_string_map():
    env = {
        "subsystems": [
            {"name": "thermal", "solver": _thermal_solver, "model": {"temperature": 300.0}},
            {"name": "solid", "solver": _solid_solver, "model": {"load": 0.0}},
        ],
        "couplings": [
            {
                "src": "thermal",
                "src_field": "temperature",
                "dst": "solid",
                "dst_field": "load",
                "map": lambda temp: temp - 299.0,
            }
        ],
        "mode": "transient",
        "nsteps": 2,
        "dt": 0.5,
        "time_update": lambda step, t, subsystems: subsystems["thermal"].model.update({"temperature": 300.0 + step}),
    }
    result = solve_study_fusion(env)
    assert result["mode"] == "transient"
    assert len(result["times"]) == 3
    assert len(result["history"]["thermal"]) == 3
