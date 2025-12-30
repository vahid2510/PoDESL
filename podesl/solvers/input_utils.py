from __future__ import annotations

from typing import Iterable, List, Mapping, Sequence

_LOAD_DOF_ALIASES = {
    "fx": "ux",
    "fy": "uy",
    "fz": "uz",
    "mx": "rx",
    "my": "ry",
    "mz": "rz",
}


def normalize_dirichlet(entries: Iterable, dof_map: Mapping[str, int]) -> List[tuple[int, str, float]]:
    normalized: List[tuple[int, str, float]] = []
    for item in entries or []:
        if isinstance(item, Mapping):
            node = int(item["node"])
            dof = str(item["dof"]).lower()
            value = float(item.get("value", 0.0))
        elif isinstance(item, Sequence) and len(item) >= 3:
            node = int(item[0])
            dof = str(item[1]).lower()
            value = float(item[2])
        else:
            continue
        if dof not in dof_map:
            raise ValueError(f"Unsupported DOF '{dof}'.")
        normalized.append((node, dof, value))
    return normalized


def normalize_point_loads(entries: Iterable, dof_order: Sequence[str]) -> List[tuple[int, List[float]]]:
    normalized: List[tuple[int, List[float]]] = []
    ncomp = len(dof_order)
    for item in entries or []:
        forces = [0.0] * ncomp
        if isinstance(item, Mapping):
            node = int(item["node"])
            dof = str(item["dof"]).lower()
            dof = _LOAD_DOF_ALIASES.get(dof, dof)
            value = float(item.get("value", 0.0))
            if dof not in dof_order:
                raise ValueError(f"Unsupported load DOF '{dof}'.")
            forces[dof_order.index(dof)] = value
        elif isinstance(item, Sequence):
            node = int(item[0])
            remainder = list(item[1:])
            if not remainder:
                continue
            if isinstance(remainder[0], str):
                idx = 0
                while idx < len(remainder):
                    token = remainder[idx]
                    if not isinstance(token, str):
                        break
                    dof = _LOAD_DOF_ALIASES.get(token.lower(), token.lower())
                    if dof not in dof_order:
                        raise ValueError(f"Unsupported load DOF '{dof}'.")
                    value = float(remainder[idx + 1]) if idx + 1 < len(remainder) else 0.0
                    forces[dof_order.index(dof)] = value
                    idx += 2
                remainder = remainder[idx:]
                if not remainder:
                    normalized.append((node, forces))
                    continue
            for idx in range(min(len(remainder), ncomp)):
                forces[idx] = float(remainder[idx])
        else:
            continue
        normalized.append((node, forces))
    return normalized


def merged_dirichlet(
    *sources: Iterable,
    dof_map: Mapping[str, int],
) -> List[List[object]]:
    """
    Merge multiple BC collections (dict or tuple style) into the legacy list format
    expected by older solvers.
    """

    combined: List = []
    for src in sources:
        if not src:
            continue
        combined.extend(src)
    normalized = normalize_dirichlet(combined, dof_map)
    return [[node, dof, value] for node, dof, value in normalized]


def merged_point_loads(
    *sources: Iterable,
    dof_order: Sequence[str],
) -> List[tuple[int, List[float]]]:
    """
    Merge loads/point_loads collections and normalize component ordering.
    """

    combined: List = []
    for src in sources:
        if not src:
            continue
        combined.extend(src)
    return normalize_point_loads(combined, dof_order)
