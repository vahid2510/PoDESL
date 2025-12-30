from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Dict, Tuple

_ALIASES = {
    "solid": "solid",
    "structure": "solid",
    "truss2d": "truss2d",
    "truss": "truss2d",
    "frame": "frame2d",
    "frame2d": "frame2d",
    "beam": "beam",
    "bar": "bar1d",
    "bar1d": "bar1d",
    "axequalsb": "axequalsb",
    "axequalb": "axequalsb",
    "equations": "equations",
    "eq": "equations",
    "linearalgebra": "linearalgebra",
    "lin-algebra": "linearalgebra",
    "thermal": "thermal",
    "heat1d": "heat1d",
    "heat": "heat1d",
    "dynamics": "dynamics",
    "mcksystem": "mcksystem",
    "mck": "mcksystem",
    "axisymmetric": "axisym",
    "axisym": "axisym",
    "static": "static",
    "steady": "steady",
    "transient": "transient",
    "direct": "direct",
}

_DOCUMENTED: Dict[Tuple[str, str, str], str] = {
    ("solid", "truss2d", "static"): "SOLID.Truss2D.Static",
    ("solid", "frame2d", "static"): "SOLID.Frame2D.Static",
    ("solid", "beam", "static"): "SOLID.Beam.Static",
    ("solid", "bar1d", "static"): "SOLID.Bar1D.Static",
    ("linearalgebra", "axequalsb", "direct"): "LA.AxEq.Direct",
    ("linearalgebra", "equations", "direct"): "LA.Equations.Direct",
    ("dynamics", "mcksystem", "transient"): "DYN.MCK.Transient",
    ("thermal", "heat1d", "steady"): "TH.Heat1D.Steady",
}


def _norm(text: str) -> str:
    if text is None:
        return ""
    s = (
        text.strip()
        .lower()
        .replace(" ", "")
        .replace("-", "")
        .replace("_", "")
        .replace("\ufeff", "")
        .replace(".", "")
    ).strip()
    return _ALIASES.get(s, s)


def _camelize(token: str) -> str:
    pieces = re.findall(r"[a-zA-Z]+|\d+", token)
    out = []
    for part in pieces:
        if part.isdigit():
            out.append(part)
        else:
            out.append(part.capitalize())
    return "".join(out)


@dataclass(frozen=True)
class ProblemDescriptor:
    title: str
    domain: str
    categories: Tuple[str, ...]
    objective: str
    raw_categories: Tuple[str, ...]
    raw_objective: str

    @property
    def primary_category(self) -> str:
        return self.categories[0] if self.categories else ""

    @property
    def category_key(self) -> str:
        return ".".join(self.categories)

    @property
    def code(self) -> str:
        domain = self.domain.upper()
        middle = ".".join(_camelize(cat) for cat in self.categories if cat)
        objective = _camelize(self.objective)
        if not middle:
            return f"{domain}.{objective}"
        return f"{domain}.{middle}.{objective}"


def parse_problem_title(title: str) -> ProblemDescriptor:
    parts = [p.strip() for p in (title or "").split(":") if p.strip()]
    if len(parts) < 3:
        raise ValueError("Invalid PROBLEM header")
    domain = _norm(parts[0])
    objective = _norm(parts[-1])
    middle_raw = [p for p in parts[1:-1] if p]
    normalized_middle = [_norm(p) for p in middle_raw if _norm(p)]
    if not normalized_middle:
        raise ValueError("Invalid PROBLEM header (missing category)")
    return ProblemDescriptor(
        title=title,
        domain=domain,
        categories=tuple(normalized_middle),
        objective=objective,
        raw_categories=tuple(middle_raw),
        raw_objective=parts[-1],
    )


def _documented_key(descriptor: ProblemDescriptor) -> Tuple[str, str, str]:
    return descriptor.domain, descriptor.category_key, descriptor.objective


def problem_code(title: str) -> str:
    descriptor = parse_problem_title(title)
    key = _documented_key(descriptor)
    if key in _DOCUMENTED:
        return _DOCUMENTED[key]
    return descriptor.code


def is_supported(title: str) -> bool:
    try:
        descriptor = parse_problem_title(title)
    except Exception:
        return False
    return _documented_key(descriptor) in _DOCUMENTED
