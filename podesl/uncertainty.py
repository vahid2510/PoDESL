from __future__ import annotations

import copy
from typing import Any, Callable, Dict, Iterable

import numpy as np

__all__ = ["sample_parameters", "mc_propagate"]


_SUPPORTED_DISTS = {"normal", "uniform", "lognormal"}


def _rng(rng):
    return rng if rng is not None else np.random.default_rng()


def sample_parameters(param_defs: Dict[str, Dict[str, Any]], rng: np.random.Generator | None = None) -> Dict[str, Any]:
    """Draw a single sample of uncertain parameters."""

    generator = _rng(rng)
    sample: Dict[str, Any] = {}
    for name, spec in (param_defs or {}).items():
        if spec is None:
            continue
        dist = str(spec.get("dist", "normal")).lower()
        if dist not in _SUPPORTED_DISTS:
            raise ValueError(f"Unsupported distribution '{dist}' for parameter '{name}'.")
        if dist == "normal":
            mean = float(spec.get("mean", 0.0))
            std = float(spec.get("std", spec.get("sigma", 1.0)))
            sample[name] = generator.normal(loc=mean, scale=std)
        elif dist == "uniform":
            low = float(spec.get("low", spec.get("min", 0.0)))
            high = float(spec.get("high", spec.get("max", 1.0)))
            sample[name] = generator.uniform(low, high)
        elif dist == "lognormal":
            mu = float(spec.get("mu", spec.get("mean", 0.0)))
            sigma = float(spec.get("sigma", spec.get("std", 1.0)))
            sample[name] = generator.lognormal(mean=mu, sigma=sigma)
    return sample


def _update_model(base: Dict[str, Any], updates: Dict[str, Any]) -> Dict[str, Any]:
    model = copy.deepcopy(base)
    for key, value in (updates or {}).items():
        model[key] = value
    return model


def _extract_path(data: Dict[str, Any], path: str) -> Any:
    if not path:
        return None
    current = data
    for token in path.split('.'):
        if isinstance(current, dict) and token in current:
            current = current[token]
        else:
            raise KeyError(f"Output '{path}' not found in solver result.")
    return current


def mc_propagate(
    solver_func: Callable[[Dict[str, Any]], Dict[str, Any]],
    base_model: Dict[str, Any],
    param_defs: Dict[str, Dict[str, Any]],
    outputs: Iterable[str],
    nsamples: int,
    rng_seed: int | None = None,
) -> Dict[str, Any]:
    """Monte Carlo uncertainty propagation driver."""

    if nsamples <= 0:
        raise ValueError("nsamples must be positive")
    outputs = list(outputs or [])
    rng = np.random.default_rng(rng_seed)

    raw_data: Dict[str, list] = {name: [] for name in outputs}

    for _ in range(nsamples):
        sample = sample_parameters(param_defs, rng)
        model = _update_model(base_model, sample)
        result = solver_func(model)
        for name in outputs:
            value = _extract_path(result, name)
            raw_data[name].append(np.asarray(value))

    stats = {"samples": nsamples, "mean": {}, "std": {}, "min": {}, "max": {}, "raw": {}}

    for name in outputs:
        values = np.stack(raw_data[name], axis=0)
        stats["raw"][name] = values
        stats["mean"][name] = values.mean(axis=0)
        stats["std"][name] = values.std(axis=0, ddof=0)
        stats["min"][name] = values.min(axis=0)
        stats["max"][name] = values.max(axis=0)

    return stats
