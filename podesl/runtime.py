from __future__ import annotations

import math
from functools import lru_cache
import keyword
from typing import Any, Callable, Dict, Iterable, Optional

import numpy as np


def report_print(*args):
    print(*args)


def report_export_csv(path: str, *arrays):
    if str(path).lower().endswith(".vtk"):
        payload = arrays[0] if arrays else None
        _export_vtk(path, payload)
        return

    cols = []
    headers = []
    for arr in arrays:
        vector = np.asarray(arr).ravel() if not isinstance(arr, dict) else None
        if vector is not None:
            cols.append(vector)
            headers.append("v" + str(len(headers) + 1))
        else:
            for key, val in arr.items():
                vv = np.asarray(val).ravel()
                cols.append(vv)
                headers.append(str(key))
    n = max((len(c) for c in cols), default=0)
    matrix = np.zeros((n, len(cols)))
    for idx, col in enumerate(cols):
        matrix[: len(col), idx] = col
    np.savetxt(path, matrix, delimiter=",", header=",".join(headers), comments="")
    print(f"exported -> {path}")


def report_plot(*series, title=None):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - matplotlib optional
        print("matplotlib not available:", exc)
        return

    i = 0
    while i < len(series):
        if i + 1 < len(series):
            x = np.ravel(series[i])
            y = np.ravel(series[i + 1])
            if x.size == y.size:
                plt.plot(x, y)
                i += 2
                continue
        plt.plot(np.ravel(series[i]))
        i += 1
    if title:
        plt.title(title)
    plt.grid(True)
    plt.show()


def _call_linspace(start, stop, num=50):
    return np.linspace(start, stop, int(num))


def _call_array(data):
    return np.array(data, dtype=float)


def _call_zeros(n):
    return np.zeros(int(n))


def _call_ones(n):
    return np.ones(int(n))


def _call_eye(n):
    return np.eye(int(n))


def _call_diag(entries: Iterable[float]):
    return np.diag(np.asarray(entries, dtype=float))


def _call_help(title_or_code: str) -> str:
    from .help import help_problem

    try:
        return help_problem(title_or_code)
    except Exception as exc:
        return f"help unavailable: {exc}"


_SAFE_CALL_TABLE: Dict[str, Callable[..., object]] = {
    "linspace": _call_linspace,
    "array": _call_array,
    "zeros": _call_zeros,
    "ones": _call_ones,
    "eye": _call_eye,
    "diag": _call_diag,
    "help": _call_help,
}


SAFE_CALLS: Dict[str, Callable[..., object]] = dict(_SAFE_CALL_TABLE)

_EVAL_GLOBALS: Dict[str, Any] = {"__builtins__": {}, "np": np, "math": math, "None": None}
_EVAL_GLOBALS.update(
    {
        "true": True,
        "false": False,
        "True": True,
        "False": False,
    }
)
_EVAL_GLOBALS.update(SAFE_CALLS)


@lru_cache(maxsize=1024)
def _compile_expr(expr: str):
    return compile(expr, "<podesl-expr>", "eval")


def eval_expr(expr: str, env: Dict[str, Any], result: Optional[Dict[str, Any]] = None):
    if expr is None:
        raise ValueError("Expression is missing.")
    text = expr.strip()
    if not text:
        raise ValueError("Expression is empty.")
    locals_dict: Dict[str, Any] = dict(env)
    if result:
        locals_dict.update(result)
        locals_dict.setdefault("result", result)
        locals_dict.setdefault("vtk", result)
    if keyword.iskeyword(text) and text in locals_dict:
        return locals_dict[text]
    code = _compile_expr(text)
    return eval(code, _EVAL_GLOBALS, locals_dict)


def lookup_report_value(name: str, env: Dict[str, Any], result: Optional[Dict[str, Any]] = None):
    if result and name in result:
        return result[name]
    if name == "vtk" and result is not None:
        return result
    if env and name in env:
        return env[name]
    raise NameError(name)


def _export_vtk(path: str, payload: Any) -> None:
    if payload is None:
        raise ValueError("VTK export requires a data object.")
    if not isinstance(payload, dict):
        raise TypeError("VTK export expects a dict containing nodes/elems.")
    nodes = payload.get("nodes")
    elems = payload.get("elems")
    if nodes is None or elems is None:
        raise KeyError("VTK payload must include 'nodes' and 'elems'.")
    nodes = np.asarray(nodes, dtype=float)
    elems = np.asarray(elems, dtype=int)
    npts = nodes.shape[0]
    ncells = elems.shape[0]
    point_data: Dict[str, np.ndarray] = {}
    cell_data: Dict[str, np.ndarray] = {}
    for key, val in payload.items():
        if key in {"nodes", "elems"}:
            continue
        _maybe_record_field(point_data, key, val, npts)
        _maybe_record_field(cell_data, key, val, ncells)
    from .solvers.vtk_writer import write_vtk_any

    write_vtk_any(path, nodes, elems, point_data or None, cell_data or None)
    print(f"exported -> {path}")


def _maybe_record_field(container: Dict[str, np.ndarray], name: str, value: Any, size: int) -> None:
    if size <= 0:
        return
    try:
        arr = np.asarray(value, dtype=float)
    except Exception:
        return
    if arr.ndim == 0 or arr.shape[0] != size:
        return
    if arr.ndim == 1:
        container[name] = arr
        return
    if arr.ndim == 2 and arr.shape[1] in (2, 3):
        container[name] = arr
        return
    if arr.ndim >= 2:
        flat = arr.reshape(size, -1)
        if flat.shape[1] in (2, 3):
            container[name] = flat
            return
        for idx in range(flat.shape[1]):
            container[f"{name}_{idx}"] = flat[:, idx]
