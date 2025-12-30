from __future__ import annotations

from typing import Dict, Iterable, Mapping, Sequence, Tuple, Union

import numpy as np

Number = Union[int, float]


def _ensure_tuple(target: Union[int, Sequence[int]]) -> Tuple[int, ...]:
    if isinstance(target, Sequence) and not isinstance(target, (str, bytes)):
        return tuple(int(idx) for idx in target)
    return (int(target),)


def _normalized_map(dof_map: Mapping[str, Union[int, Sequence[int]]]) -> Dict[str, Tuple[int, ...]]:
    return {key.lower(): _ensure_tuple(val) for key, val in dof_map.items()}


def build_dirichlet_dict(
    entries: Iterable[Sequence[object]] | None,
    *,
    ndof_per_node: int,
    dof_map: Mapping[str, Union[int, Sequence[int]]],
) -> Dict[int, float]:
    """Return a mapping of constrained global DOF indices to prescribed values."""

    mapping = _normalized_map(dof_map)
    fixed: Dict[int, float] = {}
    if not entries:
        return fixed

    for raw_entry in entries:
        if raw_entry is None:
            continue
        entry = list(raw_entry)
        if not entry:
            continue
        node = int(entry[0])
        rest = entry[1:]
        if not rest:
            continue
        first = rest[0]
        if isinstance(first, str):
            name = first.strip().lower()
            indices = mapping.get(name)
            if indices is None:
                raise ValueError(f"Unsupported DOF '{first}' in fix entry.")
            value = float(rest[1]) if len(rest) >= 2 else 0.0
            for local in indices:
                idx = ndof_per_node * node + int(local)
                fixed[idx] = value
            continue

        value = float(first)
        for local in range(ndof_per_node):
            idx = ndof_per_node * node + local
            fixed[idx] = value

    return fixed


def dirichlet_index_values(
    entries: Iterable[Sequence[object]] | None,
    *,
    ndof_per_node: int,
    dof_map: Mapping[str, Union[int, Sequence[int]]],
) -> Tuple[np.ndarray, np.ndarray]:
    fixed = build_dirichlet_dict(entries, ndof_per_node=ndof_per_node, dof_map=dof_map)
    if not fixed:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=float)
    idx = np.array(sorted(fixed.keys()), dtype=int)
    vals = np.array([fixed[i] for i in idx], dtype=float)
    return idx, vals


def dirichlet_vector(
    entries: Iterable[Sequence[object]] | None,
    *,
    ndof: int,
    ndof_per_node: int,
    dof_map: Mapping[str, Union[int, Sequence[int]]],
) -> Tuple[np.ndarray, np.ndarray]:
    idx, vals = dirichlet_index_values(
        entries, ndof_per_node=ndof_per_node, dof_map=dof_map
    )
    vec = np.zeros(ndof, dtype=float)
    if idx.size:
        vec[idx] = vals
    return idx, vec


PLANAR_UV_DOF_MAP: Dict[str, Tuple[int, ...]] = {
    "ux": (0,),
    "x": (0,),
    "u": (0,),
    "ur": (0,),
    "uy": (1,),
    "y": (1,),
    "v": (1,),
    "both": (0, 1),
    "uv": (0, 1),
    "all": (0, 1),
}

FRAME2D_DOF_MAP: Dict[str, Tuple[int, ...]] = {
    "ux": (0,),
    "x": (0,),
    "uy": (1,),
    "y": (1,),
    "rz": (2,),
    "theta": (2,),
    "rot": (2,),
    "rotation": (2,),
    "both": (0, 1, 2),
    "uv": (0, 1),
    "all": (0, 1, 2),
    "clamp": (0, 1, 2),
    "clamped": (0, 1, 2),
}

FRAME3D_DOF_MAP: Dict[str, Tuple[int, ...]] = {
    "ux": (0,),
    "x": (0,),
    "uy": (1,),
    "y": (1,),
    "uz": (2,),
    "z": (2,),
    "rx": (3,),
    "theta_x": (3,),
    "mx": (3,),
    "ry": (4,),
    "theta_y": (4,),
    "my": (4,),
    "rz": (5,),
    "theta_z": (5,),
    "mz": (5,),
    "both": (0, 1, 2),
    "all": (0, 1, 2, 3, 4, 5),
    "clamp": (0, 1, 2, 3, 4, 5),
    "clamped": (0, 1, 2, 3, 4, 5),
}

SOLID3D_DOF_MAP: Dict[str, Tuple[int, ...]] = {
    "ux": (0,),
    "x": (0,),
    "uy": (1,),
    "y": (1,),
    "uz": (2,),
    "z": (2,),
    "all": (0, 1, 2),
    "both": (0, 1, 2),
    "clamp": (0, 1, 2),
    "clamped": (0, 1, 2),
}
