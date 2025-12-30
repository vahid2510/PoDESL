import numpy as np

from podesl.solvers.bar1d import solve_bar1d_static
from podesl.solvers.truss2d import solve_truss2d
from podesl.solvers.solid3d_hex8_transient import solve_solid3d_hex8_transient
from podesl.solvers.solid3d_hex8_nonlinear import solve_solid3d_hex8_nonlinear
from podesl.solvers.solid3d_hex20 import solve_solid3d_hex20_static, solve_solid3d_hex20_modal
from podesl.solvers.frame2d_beamcolumn import solve_frame2d_beamcolumn_static, solve_frame2d_beamcolumn_nonlinear
from podesl.solvers.frame3d_beamcolumn import solve_frame3d_beamcolumn_static, solve_frame3d_beamcolumn_modal
from podesl.solvers.shell_membrane_static import solve_shell_membrane_static
from podesl.solvers.thermal1d import solve_thermal_1d_steady, solve_thermal_1d_transient
from podesl.solvers.thermomech3d_coupled import solve_thermomechanical_3d_coupled
from podesl.solvers.fsi2d_coupled import solve_fsi2d_coupled


def test_bar1d_static_sigma_shape():
    res = solve_bar1d_static(
        {
            "L": 1.0,
            "A": 1.0e-4,
            "E": 210e9,
            "nel": 4,
            "body_force": 0.0,
            "left": "fixed",
            "right": "traction",
            "traction_right": 1e5,
        }
    )
    assert res["sigma"].shape == (4,)


def test_truss2d_static_displacement_nonzero():
    res = solve_truss2d(
        {
            "nodes": np.array([[0, 0], [1, 0], [0.5, 0.8]]),
            "elems": np.array([[0, 2], [1, 2], [0, 1]]),
            "E": 210e9,
            "A": np.array([1e-4, 1e-4, 1e-4]),
            "bcs": [
                {"node": 0, "dof": "ux", "value": 0.0},
                {"node": 0, "dof": "uy", "value": 0.0},
                {"node": 1, "dof": "uy", "value": 0.0},
            ],
            "loads": [
                {"node": 2, "dof": "uy", "value": -800.0},
                [2, 0.0, -200.0],
            ],
        }
    )
    assert np.isfinite(res["U"]).all()


def test_solid3d_hex8_transient_shapes():
    nodes = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    )
    elements = np.array([[0, 1, 2, 3, 4, 5, 6, 7]])

    def _zero_load(_):
        return np.zeros(24)

    res = solve_solid3d_hex8_transient(
        {
            "nodes": nodes,
            "elements": elements,
            "E": 210e9,
            "nu": 0.3,
            "rho": 7800.0,
            "damping_ratio": 0.02,
            "dt": 1e-4,
            "t_end": 2e-4,
            "fixed_dofs": [0, 1, 2],
            "load_func": _zero_load,
        }
    )
    assert res["u"].shape[0] == 3


def test_solid3d_hex8_nonlinear_u_size():
    nodes = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    )
    elements = np.array([[0, 1, 2, 3, 4, 5, 6, 7]])
    res = solve_solid3d_hex8_nonlinear(
        {
            "nodes": nodes,
            "elements": elements,
            "material": {"E": 210e9, "nu": 0.3, "nonlinear_coef": 1e5},
            "loads": [{"node": 6, "dof": "uy", "value": -1000.0}],
            "bcs": [
                {"node": 0, "dof": "ux", "value": 0.0},
                {"node": 0, "dof": "uy", "value": 0.0},
                {"node": 0, "dof": "uz", "value": 0.0},
            ],
            "max_iter": 5,
            "load_steps": 1,
            "tol": 1e-3,
        }
    )
    assert res["u"].shape == (24,)


def _hex20_mesh():
    nodes = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [0.5, 0, 0],
            [1, 0.5, 0],
            [0.5, 1, 0],
            [0, 0.5, 0],
            [0, 0, 0.5],
            [1, 0, 0.5],
            [1, 1, 0.5],
            [0, 1, 0.5],
            [0.5, 0, 1],
            [1, 0.5, 1],
            [0.5, 1, 1],
            [0, 0.5, 1],
        ]
    )
    elements = np.arange(20, dtype=int).reshape(1, 20)
    return nodes, elements


def test_solid3d_hex20_static_runs():
    nodes, elements = _hex20_mesh()
    res = solve_solid3d_hex20_static(
        {
            "nodes": nodes,
            "elements": elements,
            "E": 210e9,
            "nu": 0.3,
            "loads": [{"node": 6, "dof": "uz", "value": -1000.0}],
            "bcs": [{"node": 0, "dof": "ux", "value": 0.0}],
        }
    )
    assert res["u"].shape[0] == 60


def test_solid3d_hex20_modal_returns_freq():
    nodes, elements = _hex20_mesh()
    res = solve_solid3d_hex20_modal(
        {
            "nodes": nodes,
            "elements": elements,
            "E": 210e9,
            "nu": 0.3,
            "rho": 7800.0,
            "nmodes": 2,
            "bcs": [{"node": 0, "dof": "ux", "value": 0.0}],
        }
    )
    assert res["freq"].shape[0] == 2


def test_shell_membrane_static_runs():
    model = {
        "nodes": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        "elements": np.array([[0, 1, 2]]),
        "thickness": 0.01,
        "E": 210e9,
        "nu": 0.3,
        "bcs": [
            {"node": 0, "dof": "ux", "value": 0.0},
            {"node": 0, "dof": "uy", "value": 0.0},
            {"node": 1, "dof": "uy", "value": 0.0},
        ],
        "loads": [
            {"node": 2, "dof": "uy", "value": -1e3},
        ],
    }
    res = solve_shell_membrane_static(model)
    assert res["u"].shape[0] == 6
    assert np.isfinite(res["u"]).all()


def test_frame2d_beamcolumn_static_disp():
    res = solve_frame2d_beamcolumn_static(
        {
            "nodes": np.array([[0, 0], [0, 3]]),
            "elements": np.array([[0, 1]]),
            "E": 210e9,
            "A": 4e-4,
            "I": 8e-6,
            "fix": [[0, "ux", 0.0], [0, "uy", 0.0], [0, "rz", 0.0]],
            "loads": [[1, 0.0, -1000.0, 0.0]],
        }
    )
    assert np.isfinite(res["disp"]).all()


def test_frame2d_beamcolumn_nonlinear_runs():
    res = solve_frame2d_beamcolumn_nonlinear(
        {
            "nodes": np.array([[0, 0], [0, 3]]),
            "elements": np.array([[0, 1]]),
            "E": 210e9,
            "A": 4e-4,
            "I": 8e-6,
            "fix": [[0, "ux", 0.0], [0, "uy", 0.0], [0, "rz", 0.0]],
            "loads": [[1, 0.0, -500.0, 0.0]],
            "load_steps": 2,
            "max_iter": 5,
        }
    )
    assert np.isfinite(res["U"]).all()


def test_frame3d_beamcolumn_static_runs():
    res = solve_frame3d_beamcolumn_static(
        {
            "nodes": np.array([[0, 0, 0], [0, 0, 3]]),
            "elements": np.array([[0, 1]]),
            "E": 210e9,
            "G": 80e9,
            "A": 4e-4,
            "Iy": 8e-6,
            "Iz": 8e-6,
            "J": 1e-5,
            "bcs": [{"node": 0, "dof": dof, "value": 0.0} for dof in ["ux", "uy", "uz", "rx", "ry", "rz"]],
            "loads": [{"node": 1, "dof": "uy", "value": -500.0}],
            "point_loads": [[1, 0.0, -500.0, 0.0, 0.0, 0.0]],
        }
    )
    assert res["u"].shape[0] == 12


def test_frame3d_beamcolumn_modal_freq():
    res = solve_frame3d_beamcolumn_modal(
        {
            "nodes": np.array([[0, 0, 0], [0, 0, 3]]),
            "elements": np.array([[0, 1]]),
            "E": 210e9,
            "G": 80e9,
            "A": 4e-4,
            "Iy": 8e-6,
            "Iz": 8e-6,
            "J": 1e-5,
            "rho": 7800.0,
            "bcs": [{"node": 0, "dof": dof, "value": 0.0} for dof in ["ux", "uy", "uz", "rx", "ry", "rz"]],
            "nmodes": 2,
        }
    )
    assert res["freq"].shape[0] == 2


def test_thermal1d_steady_profile():
    res = solve_thermal_1d_steady({"L": 1.0, "A": 1.0, "k": 10.0, "nel": 4, "T_left": 100.0, "T_right": 50.0})
    assert np.isfinite(res["T"]).all()


def test_thermal1d_transient_history():
    res = solve_thermal_1d_transient(
        {
            "L": 1.0,
            "A": 1.0,
            "k": 10.0,
            "rho": 7800.0,
            "c": 500.0,
            "nel": 4,
            "dt": 0.1,
            "t_end": 0.2,
            "T_left": 100.0,
            "T_right": 50.0,
        }
    )
    assert res["T"].shape[0] == 3


def test_thermomech3d_coupled_returns_histories():
    nodes = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ]
    )
    elements = np.array([[0, 1, 2, 3, 4, 5, 6, 7]])
    thermal_model = {
        "nodes": nodes,
        "elems": elements,
        "k": 10.0,
        "rho": 7800.0,
        "c": 500.0,
        "dt": 0.1,
        "t_end": 0.2,
        "store_every": 1,
        "fixT": [[0, 100.0]],
        "T0": 20.0,
    }
    mechanical_model = {
        "nodes": nodes,
        "elements": elements,
        "E": 210e9,
        "nu": 0.3,
        "loads": [],
        "bcs": [{"node": 0, "dof": "ux", "value": 0.0}],
    }
    res = solve_thermomechanical_3d_coupled(
        {
            "thermal_model": thermal_model,
            "mechanical_model": mechanical_model,
            "alpha": 1e-5,
            "T_ref": 20.0,
            "dt": 0.1,
            "t_end": 0.2,
        }
    )
    assert res["displacement"].shape[0] == 3


def test_fsi2d_coupled_history():
    structure = {
        "nodes": np.array([[0, 0], [0, 3]]),
        "elems": np.array([[0, 1]]),
        "E": 210e9,
        "A": 4e-4,
        "I": 8e-6,
        "fix": [[0, "ux", 0.0], [0, "uy", 0.0], [0, "rz", 0.0]],
        "pressure_nodes": [1],
    }

    def pressure(t):
        return 1000.0 * np.sin(2 * np.pi * t)

    res = solve_fsi2d_coupled({"structure_model": structure, "pressure_time_func": pressure, "dt": 0.1, "t_end": 0.2})
    assert res["disp"].shape[0] == 3
