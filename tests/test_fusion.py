import numpy as np

from podesl.fusion import CouplingLink, FusionEngine, Subsystem


def _thermal_solver(model):
    temp = model.get("temperature", 300.0)
    return {"temperature": temp}


def _solid_solver(model):
    load = model.get("loads", {}).get("thermal", 0.0)
    return {"tip_disp": load * 1e-3}


def _first_subsystem():
    return Subsystem("thermal", _thermal_solver, {"temperature": 350.0})


def _second_subsystem():
    return Subsystem("solid", _solid_solver, {"loads": {"thermal": 0.0}})


def test_subsystem_field_access():
    thermal = _first_subsystem()
    thermal.solve()
    assert thermal.get_field("temperature") == 350.0
    solid = _second_subsystem()
    solid.set_model_field("loads.thermal", 42.0)
    assert solid.model["loads"]["thermal"] == 42.0


def test_fusion_run_steady():
    thermal = _first_subsystem()
    solid = _second_subsystem()

    link = CouplingLink(
        src="thermal",
        src_field="temperature",
        dst="solid",
        dst_field="loads.thermal",
        map_func=lambda temp: 1.2e-5 * (temp - 300.0),
    )

    engine = FusionEngine([thermal, solid], [link], max_iter=3, tol=1e-6)
    results = engine.run_steady()
    assert "thermal" in results and "solid" in results
    assert solid.model["loads"]["thermal"] == 1.2e-5 * 50.0
    assert np.isclose(results["solid"]["tip_disp"], solid.model["loads"]["thermal"] * 1e-3)


def test_fusion_run_transient():
    calls = []

    def time_update(step, t, subsystems):
        subsystems["thermal"].model["temperature"] = 300.0 + 10.0 * step
        calls.append((step, t))

    thermal = _first_subsystem()
    solid = _second_subsystem()

    link = CouplingLink(
        src="thermal",
        src_field="temperature",
        dst="solid",
        dst_field="loads.thermal",
        map_func=lambda temp: temp - 250.0,
    )

    engine = FusionEngine([thermal, solid], [link])
    times, history = engine.run_transient(nsteps=3, time_update=time_update, dt=0.1)
    assert len(times) == 4
    assert len(history["thermal"]) == 4
    assert calls[0][0] == 0
