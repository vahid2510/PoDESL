from __future__ import annotations

from typing import Callable, Dict, List

import numpy as np

__all__ = ["finite_diff_gradient", "optimize_design"]


def finite_diff_gradient(f: Callable[[np.ndarray], float], x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    grad = np.zeros_like(x)
    for i in range(x.size):
        step = eps * max(1.0, abs(x[i]))
        x_forward = x.copy()
        x_backward = x.copy()
        x_forward[i] += step
        x_backward[i] -= step
        f_forward = float(f(x_forward))
        f_backward = float(f(x_backward))
        grad[i] = (f_forward - f_backward) / (2.0 * step)
    return grad


def optimize_design(
    f_obj: Callable[[np.ndarray], float],
    x0: np.ndarray,
    max_iter: int = 100,
    step_size: float = 1e-2,
    tol: float = 1e-6,
    verbose: bool = False,
) -> Dict[str, object]:
    x = np.asarray(x0, dtype=float).copy()
    history: Dict[str, List] = {"x": [], "f": [], "grad_norm": []}
    converged = False

    f_current = float(f_obj(x))
    history["x"].append(x.copy())
    history["f"].append(f_current)
    history["grad_norm"].append(np.nan)

    iterations = 0
    for iteration in range(1, max_iter + 1):
        iterations = iteration
        grad = finite_diff_gradient(f_obj, x)
        grad_norm = float(np.linalg.norm(grad))
        history["x"].append(x.copy())
        history["f"].append(f_current)
        history["grad_norm"].append(grad_norm)

        if verbose:
            print(f"iter {iteration}: f={f_current:.6e}, |grad|={grad_norm:.3e}")
        if grad_norm < tol:
            converged = True
            break

        x_new = x - step_size * grad
        f_new = float(f_obj(x_new))

        if not np.isfinite(f_new) or abs(f_new) > 1e20:
            break

        x = x_new
        f_current = f_new

    return {
        "x_opt": x,
        "f_opt": f_current,
        "history": history,
        "converged": converged,
        "iterations": iterations,
    }
