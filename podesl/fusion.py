from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

__all__ = ["Subsystem", "CouplingLink", "FusionEngine"]


def _split_path(path: str) -> List[str]:
    if not path:
        raise ValueError("Field path must be non-empty")
    return [segment.strip() for segment in path.split(".") if segment.strip()]


def _get_from_path(container: Dict[str, Any], path: str) -> Any:
    current: Any = container
    for key in _split_path(path):
        if not isinstance(current, dict) or key not in current:
            raise KeyError(f"Missing key '{key}' while accessing '{path}'")
        current = current[key]
    return current


def _set_in_path(container: Dict[str, Any], path: str, value: Any) -> None:
    keys = _split_path(path)
    current = container
    for key in keys[:-1]:
        if key not in current or not isinstance(current[key], dict):
            current[key] = {}
        current = current[key]
    current[keys[-1]] = value


def _value_norm(value: Any) -> float:
    try:
        arr = np.asarray(value, dtype=float)
        return float(np.linalg.norm(arr))
    except (TypeError, ValueError):
        return float(value) if isinstance(value, (int, float)) else 0.0


def _diff_norm(old: Any, new: Any) -> float:
    if old is None:
        return _value_norm(new)
    try:
        arr_old = np.asarray(old, dtype=float)
        arr_new = np.asarray(new, dtype=float)
        return float(np.linalg.norm(arr_new - arr_old))
    except (TypeError, ValueError):
        return 0.0 if old == new else 1.0


class Subsystem:
    """Represents a single-physics solver in the fusion framework.

    Example
    -------
    >>> solid = Subsystem("solid", solve_solid, {"loads": {}})
    >>> solid.set_model_field("loads.thermal", 10.0)
    >>> result = solid.solve()
    >>> disp_tip = solid.get_field("results.tip_disp")
    """

    def __init__(
        self,
        name: str,
        solve_func: Callable[[Dict[str, Any]], Dict[str, Any]],
        model: Dict[str, Any],
    ) -> None:
        self.name = name
        self.solve_func = solve_func
        self.model: Dict[str, Any] = copy.deepcopy(model)
        self.last_result: Optional[Dict[str, Any]] = None

    def solve(self) -> Dict[str, Any]:
        """Run the subsystem solver and store the latest result."""

        self.last_result = self.solve_func(self.model)
        return self.last_result

    def get_field(self, field_path: str) -> Any:
        """Return a field from the last result using dotted path notation."""

        if self.last_result is None:
            raise RuntimeError(f"Subsystem '{self.name}' has not been solved yet.")
        try:
            return _get_from_path(self.last_result, field_path)
        except KeyError as exc:
            raise RuntimeError(
                f"Field '{field_path}' not found in results of subsystem '{self.name}'."
            ) from exc

    def set_model_field(self, field_path: str, value: Any) -> None:
        """Set (and create if necessary) a model entry via dotted path notation."""

        _set_in_path(self.model, field_path, value)


@dataclass
class CouplingLink:
    """Defines a unidirectional coupling between two subsystems."""

    src: str
    src_field: str
    dst: str
    dst_field: str
    map_func: Callable[[Any], Any]
    description: str = ""


class FusionEngine:
    """Couples several subsystems using sequential or fixed-point iterations.

    Example
    -------
    >>> thermal = Subsystem("thermal", solve_thermal, {"boundary": {"flux": 0.0}})
    >>> solid = Subsystem("solid", solve_solid, {"loads": {}})
    >>> link = CouplingLink(
    ...     src="thermal",
    ...     src_field="boundary.temperature",
    ...     dst="solid",
    ...     dst_field="loads.thermal_strain",
    ...     map_func=lambda temp: 1.2e-5 * (temp - 293.15),
    ... )
    >>> engine = FusionEngine([thermal, solid], [link], mode="steady")
    >>> results = engine.run_steady()
    """

    def __init__(
        self,
        subsystems: List[Subsystem],
        couplings: List[CouplingLink],
        mode: str = "steady",
        max_iter: int = 10,
        tol: float = 1e-3,
        verbose: bool = False,
    ) -> None:
        self.mode = mode
        self.max_iter = max_iter
        self.tol = tol
        self.verbose = verbose
        self._subsystems: Dict[str, Subsystem] = {}
        for subsystem in subsystems:
            if subsystem.name in self._subsystems:
                raise ValueError(f"Duplicate subsystem name '{subsystem.name}'")
            self._subsystems[subsystem.name] = subsystem
        self.couplings = couplings
        self._last_residual: float = 0.0
        self.residual_history: List[float] = []

    def get_subsystem(self, name: str) -> Subsystem:
        """Return the subsystem by name or raise KeyError."""

        return self._subsystems[name]

    def _exchange_coupled_fields(self) -> None:
        residual = 0.0
        for link in self.couplings:
            src = self.get_subsystem(link.src)
            dst = self.get_subsystem(link.dst)
            if src.last_result is None:
                raise RuntimeError(
                    f"Source subsystem '{link.src}' has no results for coupling."
                )
            value = src.get_field(link.src_field)
            mapped_value = link.map_func(value)
            old_value = None
            try:
                old_value = _get_from_path(dst.model, link.dst_field)
            except KeyError:
                old_value = None
            residual += _diff_norm(old_value, mapped_value)
            dst.set_model_field(link.dst_field, mapped_value)
        self._last_residual = residual

    def _compute_interface_residual(self) -> float:
        return self._last_residual

    def run_steady(self) -> Dict[str, Dict[str, Any]]:
        """Run a steady coupled simulation using fixed-point iterations."""

        self.residual_history = []
        for iteration in range(self.max_iter):
            for subsystem in self._subsystems.values():
                subsystem.solve()
            self._exchange_coupled_fields()
            residual = self._compute_interface_residual()
            self.residual_history.append(residual)
            if self.verbose:
                print(f"[Fusion] iter={iteration} residual={residual:.3e}")
            if residual < self.tol:
                break
        return {name: sub.last_result for name, sub in self._subsystems.items()}

    def run_transient(
        self,
        nsteps: int,
        time_update: Callable[[int, float, Dict[str, Subsystem]], None],
        dt: float,
    ) -> Tuple[np.ndarray, Dict[str, List[Dict[str, Any]]]]:
        """Run a simple staggered transient analysis.

        Parameters
        ----------
        nsteps : int
            Number of time steps (exclusive) in the simulation. Includes step zero.
        time_update : Callable
            Callback invoked before each step is solved; allows user modifications.
        dt : float
            Time increment between steps.
        """

        times = np.linspace(0.0, dt * nsteps, nsteps + 1)
        histories: Dict[str, List[Dict[str, Any]]] = {
            name: [] for name in self._subsystems
        }

        self.residual_history = []
        for step, t in enumerate(times):
            time_update(step, t, self._subsystems)
            for subsystem in self._subsystems.values():
                subsystem.model.setdefault("time", t)
                subsystem.solve()
            self._exchange_coupled_fields()
            self.residual_history.append(self._last_residual)
            for name, subsystem in self._subsystems.items():
                histories[name].append(copy.deepcopy(subsystem.last_result or {}))
        return times, histories
