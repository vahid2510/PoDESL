from __future__ import annotations

from typing import Any, Dict


def _beam_static_example() -> str:
    return """PROBLEM "Solid : Beam : Static"
GIVEN
  L = 2.0
  E = 210e9
  I = 8.0e-6
  nel = 40
  q = 1000.0
  left = "clamped"
  right = "clamped"
REPORT
  print w
  print reactions
"""


def _bar1d_example() -> str:
    return """PROBLEM "Solid : Bar1D : Static"
GIVEN
  L = 2.0
  E = 210e9
  A = 2.0e-4
  nel = 20
  p = 1000.0
  left = "fixed"
  right = "fixed"
REPORT
  print u
  export "bar_u.csv" x u
"""


def _truss2d_example() -> str:
    return """PROBLEM "Solid : Truss2D : Static"
GIVEN
  nodes = [[0,0],[1,0],[1,1]]
  elems = [[0,1, 1e-4, 210e9],
           [1,2, 1e-4, 210e9],
           [0,2, 1e-4, 210e9]]
  supports = [[0, 1,1], [1, 0,1]]  # [node, fix_ux, fix_uy]
  loads = [[2, 0, -1000]]
REPORT
  print U
  print reactions
"""


def _frame2d_example() -> str:
    return """PROBLEM "Solid : Frame2D : Static"
GIVEN
  nodes = [[0,0],[4,0],[4,3]]
  elems = [[0,1, 3.2e-3, 210e9, 8.0e-6, 0.0],
           [1,2, 3.2e-3, 210e9, 8.0e-6, 0.0]]
  supports = [[0,1,1,1],[2,0,1,0]]
  loads = [[1, 0.0, -15000.0, 0.0]]
REPORT
  print U
  print member_forces
"""


def _mck_example() -> str:
    return """PROBLEM "Dynamics : MCKSystem : Transient"
GIVEN
  M = [[2,0],[0,1]]
  C = [[0.02,0],[0,0.02]]
  K = [[2000,-1000],[-1000,1000]]
  u0 = [0,0]
  v0 = [0,0]
  dt = 0.001
  t_end = 0.5
  f = [0.0, 1.0]
REPORT
  print U
"""


def _axe_example() -> str:
    return """PROBLEM "LinearAlgebra : AxEqualsb : Direct"
GIVEN
  A = [[3,-1,0],[-1,2,-1],[0,-1,3]]
  b = [2,-2,2]
SOLVE
  method = "LU"
REPORT
  print x
"""


def _equations_example() -> str:
    return """PROBLEM "LinearAlgebra : Equations : Direct"
EQUATIONS
  unknowns = ["u0","u1","u2"]
  eq " 2*u0 - 1*u1      =  5"
  eq "-1*u0 + 2*u1 - 1*u2 = -1"
  eq "        -1*u1 + 2*u2 =  1"
REPORT
  print x
"""


def _heat1d_example() -> str:
    return """PROBLEM "Thermal : Heat1D : Steady"
GIVEN
  L = 1.0
  k = 200.0
  A = 1.0
  qdot = 0.0
  nel = 10
  left = "Dirichlet"; T_left = 100.0
  right = "Robin"; h_right = 10.0; Tinf_right = 20.0
REPORT
  print T
"""


HELP_INDEX: Dict[str, Dict[str, Any]] = {
    "SOLID.Beam.Static": {
        "title": "Solid : Beam : Static",
        "summary": "Euler-Bernoulli 1D beam under static loads.",
        "inputs_required": ["L", "E", "I"],
        "inputs_optional": ["nel", "q", "point_loads", "left", "right"],
        "outputs": ["w", "theta", "M", "reactions", "x"],
        "example": _beam_static_example(),
    },
    "SOLID.Bar1D.Static": {
        "title": "Solid : Bar1D : Static",
        "summary": "1D axial bar with distributed or nodal loads.",
        "inputs_required": ["L", "E", "A"],
        "inputs_optional": ["nel", "p", "point_loads", "left", "right", "U_left", "U_right"],
        "outputs": ["u", "eps", "sigma", "x", "xm"],
        "example": _bar1d_example(),
    },
    "SOLID.Truss2D.Static": {
        "title": "Solid : Truss2D : Static",
        "summary": "2D pin-jointed truss, linear small-displacement.",
        "inputs_required": ["nodes", "elems", "supports"],
        "inputs_optional": ["loads"],
        "outputs": ["U", "reactions", "member_forces"],
        "example": _truss2d_example(),
    },
    "SOLID.Frame2D.Static": {
        "title": "Solid : Frame2D : Static",
        "summary": "Planar beam-column frame with optional distributed loads.",
        "inputs_required": ["nodes", "elems", "supports"],
        "inputs_optional": ["loads", "nonlinear", "max_iter", "tol"],
        "outputs": ["U", "reactions", "member_forces"],
        "example": _frame2d_example(),
    },
    "DYN.MCK.Transient": {
        "title": "Dynamics : MCKSystem : Transient",
        "summary": "Second-order M u'' + C u' + K u = f(t) via Newmark-beta.",
        "inputs_required": ["M", "K", "dt", "t_end"],
        "inputs_optional": ["C", "u0", "v0", "f"],
        "outputs": ["U", "V", "A", "t"],
        "example": _mck_example(),
    },
    "LA.AxEq.Direct": {
        "title": "LinearAlgebra : AxEqualsb : Direct",
        "summary": "Solve A x = b using LU/Cholesky/QR fallback.",
        "inputs_required": ["A", "b"],
        "inputs_optional": ["method"],
        "outputs": ["x"],
        "example": _axe_example(),
    },
    "LA.Equations.Direct": {
        "title": "LinearAlgebra : Equations : Direct",
        "summary": "Solve linear equations given as algebraic strings.",
        "inputs_required": ["unknowns", "equations"],
        "inputs_optional": [],
        "outputs": ["x"],
        "example": _equations_example(),
    },
    "TH.Heat1D.Steady": {
        "title": "Thermal : Heat1D : Steady",
        "summary": "1D steady conduction with Dirichlet/Robin boundary conditions.",
        "inputs_required": ["L", "k"],
        "inputs_optional": ["A", "nel", "qdot", "left", "right", "T_left", "T_right", "h_left", "Tinf_left", "h_right", "Tinf_right"],
        "outputs": ["T", "x"],
        "example": _heat1d_example(),
    },
}
