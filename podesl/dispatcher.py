from __future__ import annotations

import re
from importlib import import_module
from pathlib import Path
from typing import Any, Callable, Dict, List, Sequence, Tuple, Optional

import numpy as np

from .problem_codes import parse_problem_title, problem_code


SolverFunc = Callable[[Dict[str, Any]], Dict[str, Any]]


def _entry(func: SolverFunc, required: List[str]) -> Dict[str, Any]:
    return {"func": func, "required": required}


_STATIC_SOLVER_SPECS: Dict[str, Tuple[str, str, List[str]]] = {
    "SOLID.Truss2D.Static": ("truss2d", "solve_truss2d", ["nodes", "elems"]),
    "SOLID.Frame2D.Static": ("frame2d", "solve_frame2d_static", ["nodes", "elems"]),
    "SOLID.Frame2D.Linear": ("frame2d", "solve_frame2d_static", ["nodes", "elems"]),
    "SOLID.Truss3D.Linear": ("truss3d", "solve_truss3d_static", ["nodes", "elems"]),
    "SOLID.Beam.Static": ("beam1d", "solve_beam_static", ["L", "E", "I"]),
    "SOLID.Bar1D.Static": ("bar1d", "solve_bar1d_static", ["L", "E", "A"]),
    "LA.AxEq.Direct": ("linalg", "solve_axequalsb_direct", ["A", "b"]),
    "LA.Equations.Direct": ("linalg", "solve_equations_direct", ["unknowns", "equations"]),
    "DYN.MCK.Transient": ("dynamics", "solve_mck_newmark", ["M", "K", "dt", "t_end"]),
    "TH.Heat1D.Steady": ("heat1d", "solve_heat1d_steady", ["L", "k"]),
}

_SOLVER_DIR = Path(__file__).resolve().with_name("solvers")
_SOLVER_MODULES = sorted(
    p.stem.lower() for p in _SOLVER_DIR.glob("*.py") if p.stem != "__init__"
)

_SOLVER_REGISTRY: Dict[str, Dict[str, Any]] = {}


_MANUAL_SOLVER_HINTS: Dict[str, Tuple[str, Optional[str], Optional[List[str]]]] = {
    "CONTACT.1D.Penalty": ("contact1d_penalty", "solve_contact1d_penalty", ["L", "E", "A", "P"]),
    "HEAT1D.Axisym.Linear.Transient": ("heat_axi_transient", "solve_heat_axi_transient", None),
    "HEAT1D.Conduction1D.Linear.Transient": ("heat1d_transient", "solve_heat1d_transient", None),
    "HEAT1D.Conduction2D.Steady": ("heat2d", "solve_heat2d", None),
    "HEAT1D.Conduction2D.Linear.Steady": ("heat2d", "solve_heat2d", None),
    "HEAT1D.Conduction2D.Linear.Transient": ("heat2d_transient", "solve_heat2d_transient", None),
    "HEAT1D.Conduction2D.Transient": ("heat2d_transient", "solve_heat2d_transient", None),
    "SOLID.Heat2D.Transient": ("heat2d_transient_bc", "solve_heat2d_transient_bc", None),
    "HEAT1D.Conduction3D.Transient": ("heat3d_transient", "solve_heat3d_transient", None),
    "HEAT1D.Oned.Linear.Transient": ("heat1d_transient", "solve_heat1d_transient", None),
    "HEAT1D.Planar.Linear.Transient": ("heat2d_planar_transient", "solve_heat2d_planar_transient", None),
    "HEAT1D.Threed.H8.Transient": ("heat3d_h8_transient", "solve_heat3d_h8_transient", None),
    "HEAT1D.Twod.Q4.Static": ("heat2d_q4_static", "solve_heat2d_q4_static", None),
    "MATERIALS.Vonmises.Planestressfe": ("ps_vonmises_fe", "solve_ps_vonmises_fe", None),
    "MATERIALS.Vonmises.Uniaxial": ("materials", "uniaxial_von_mises_curve", None),
    "SOLID.Frame2D.Buckling": ("frame2d_buckling", "solve_frame2d_buckling", ["nodes", "elems"]),
    "SOLID.Frame2D.Staticnl": ("frame2d_static_nl", "solve_frame2d_static_nl", None),
    "SOLID.Frame2D.Staticnonlinear.Pdelta": ("frame2d_pdelta", "solve_frame2d_pdelta", None),
    "SOLID.Frame2D.Transient": ("frame2d_transient_newmark", "solve_frame2d_transient_newmark", None),
    "SOLID.Axisym.Q4.Static": ("axisym_q4_static", "solve_axisym_q4_static", None),
    "SOLID.Bar1D.Thermoelastic.Static": ("bar1d_thermoelastic_static", "solve_bar1d_thermoelastic_static", ["nodes", "elems"]),
    "SOLID.Beam.Static": ("beam1d", "solve_beam_static", ["L", "E", "I"]),
    "SOLID.Planestrain.Linearfe": ("ps_linear_fe", "solve_ps_linear_fe", None),
    "SOLID.Planestrain.Q4.Static": ("plane_strain_q4_static", "solve_plane_strain_q4_static", None),
    "SOLID.Planestrain2D.Linear": ("plane_strain2d", "solve_plane_strain2d", None),
    "SOLID.Planestress.Q4.Static": ("plane_stress_q4_static", "solve_plane_stress_q4_static", None),
    "SOLID.Planestress2D.Linear": ("planestress2d", "solve_planestress2d", None),
    "SOLID.Platemindlin.Static": ("plate_mindlin_static", "solve_plate_mindlin_static", None),
    "SOLID.Platemindlin2D.Linear": ("plate_mindlin2d", "solve_plate_mindlin_q4", None),
    "SOLID.Shell.Flatdkt": ("shell_flat_dkt_proxy", "solve_shell_flat_dkt_static", None),
    "SOLID.Shell2D.Static": ("shell2d_mindlin_static", "solve_shell2d_mindlin_static", None),
    "SHELL.Membrane.Static": ("shell_membrane_static", "solve_shell_membrane_static", ["nodes", "elements"]),
    "STUDY.Montecarlo.Run": ("study_drivers", "solve_study_montecarlo", None),
    "STUDY.Optimize.Gradient": ("study_drivers", "solve_study_optimize", None),
    "STUDY.Scenario.Compare": ("study_drivers", "solve_study_scenarios", None),
    "STUDY.Fusion.Steady": ("study_fusion", "solve_study_fusion", None),
    "STUDY.Fusion.Transient": ("study_fusion", "solve_study_fusion", None),
    "SOLID.Truss2D.Bucklinglinear": ("truss2d_buckling_linear", "solve_truss2d_buckling_linear", None),
    "SOLID.Truss2D.Elastoplastic": ("truss2d_ep_bilin", "solve_truss2d_ep_bilin", None),
    "SOLID.Truss2D.Geomnl": ("truss2d_geom_nl", "solve_truss2d_geom_nl", None),
    "SOLID.Truss2D.Modal": ("truss2d_modal", "solve_truss2d_modal", ["nodes", "elems"]),
    "SOLID.Truss2D.Nonlinear": ("truss2d_nonlinear", "solve_truss2d_nonlinear", None),
    "SOLID.Truss3D.Modal": ("truss3d", "solve_truss3d_modal", None),
    "SOLID3D.Hex8.Transient": ("solid3d_hex8_transient", "solve_solid3d_hex8_transient", None),
    "SOLID3D.Hex8.Nonlinear": ("solid3d_hex8_nonlinear", "solve_solid3d_hex8_nonlinear", None),
    "SOLID3D.Hex20.Static": ("solid3d_hex20", "solve_solid3d_hex20_static", None),
    "SOLID3D.Hex20.Modal": ("solid3d_hex20", "solve_solid3d_hex20_modal", None),
    "FRAME2D.Beamcolumn.Static": ("frame2d_beamcolumn", "solve_frame2d_beamcolumn_static", None),
    "FRAME2D.Beamcolumn.Nonlinear": ("frame2d_beamcolumn", "solve_frame2d_beamcolumn_nonlinear", None),
    "FRAME3D.Beamcolumn.Static": ("frame3d_beamcolumn", "solve_frame3d_beamcolumn_static", None),
    "FRAME3D.Beamcolumn.Modal": ("frame3d_beamcolumn", "solve_frame3d_beamcolumn_modal", None),
    "THERMAL.1DConduction.Steady": ("thermal1d", "solve_thermal_1d_steady", None),
    "THERMAL.1DConduction.Transient": ("thermal1d", "solve_thermal_1d_transient", None),
    "THERMOMECHANICAL.3D.Coupled": ("thermomech3d_coupled", "solve_thermomechanical_3d_coupled", None),
    "FLUIDSTRUCTUREINTERACTION.2D.Coupled": ("fsi2d_coupled", "solve_fsi2d_coupled", None),
}

_SOLVER_HINTS_CACHE: Optional[Dict[str, Tuple[str, Optional[str], Optional[List[str]]]]] = None


def _get_solver_hints() -> Dict[str, Tuple[str, Optional[str], Optional[List[str]]]]:
    global _SOLVER_HINTS_CACHE
    if _SOLVER_HINTS_CACHE is not None:
        return _SOLVER_HINTS_CACHE
    hints: Dict[str, Tuple[str, Optional[str], Optional[List[str]]]] = dict(
        _MANUAL_SOLVER_HINTS
    )
    hints.update(_build_example_hints())
    _SOLVER_HINTS_CACHE = hints
    return hints


def _build_example_hints() -> Dict[str, Tuple[str, Optional[str], Optional[List[str]]]]:
    hints: Dict[str, Tuple[str, Optional[str], Optional[List[str]]]] = {}
    example_dir = Path(__file__).resolve().parents[1] / "examples"
    if not example_dir.exists():
        return hints
    for path in example_dir.glob("*.dsl"):
        code = _code_from_example(path)
        if not code or code in hints or code in _MANUAL_SOLVER_HINTS:
            continue
        module = _guess_module_from_stem(path.stem.lower())
        if module:
            hints[code] = (module, None, None)
    return hints


def _code_from_example(path: Path) -> Optional[str]:
    try:
        text = path.read_text(encoding="utf-8")
    except Exception:
        return None
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue
        if stripped.startswith("PROBLEM"):
            rest = stripped[len("PROBLEM") :].strip()
            if rest and rest[0] in ("'", '"') and rest[-1] == rest[0]:
                title = rest[1:-1]
            else:
                title = rest
            try:
                return problem_code(title)
            except Exception:
                return None
    return None


def _guess_module_from_stem(stem: str) -> Optional[str]:
    stem = stem.lower()
    best: Optional[str] = None
    for mod in _SOLVER_MODULES:
        if (
            stem == mod
            or stem.startswith(mod + "_")
            or mod.startswith(stem + "_")
            or stem.startswith(mod)
            or mod.startswith(stem)
        ):
            if best is None or len(mod) > len(best):
                best = mod
    return best


def _problem_metadata(problem: Dict[str, Any]) -> Dict[str, Any]:
    problem = problem or {}
    title = problem.get("title") or ""
    code = problem.get("code")
    categories = tuple(problem.get("categories") or [])
    objective = problem.get("objective")
    raw_categories = tuple(problem.get("raw_categories") or [])
    raw_objective = problem.get("raw_objective")
    if code and categories and objective:
        return {
            "code": code,
            "categories": categories,
            "objective": objective,
            "raw_categories": raw_categories,
            "raw_objective": raw_objective or "",
        }
    descriptor = parse_problem_title(title)
    return {
        "code": code or descriptor.code,
        "categories": categories or descriptor.categories,
        "objective": objective or descriptor.objective,
        "title": title,
        "raw_categories": raw_categories or descriptor.raw_categories,
        "raw_objective": raw_objective or descriptor.raw_objective,
    }


def _normalize_environment(env: Dict[str, Any], meta: Dict[str, Any]) -> Dict[str, Any]:
    env = dict(env or {})

    if "fix" not in env and "fixT" in env:
        env["fix"] = [list(entry) for entry in env["fixT"]]
    supports_original = env.get("supports")
    if supports_original is not None:
        env["supports"] = [list(entry) for entry in supports_original]
    if "fix" not in env and supports_original is not None:
        env["fix"] = [list(entry) for entry in supports_original]
    elif "fix" in env:
        env["fix"] = [list(entry) for entry in env["fix"]]
    if "cp" not in env and "c" in env:
        env["cp"] = env["c"]
    if "c" not in env and "cp" in env:
        env["c"] = env["cp"]
    if "rho" not in env and "density" in env:
        env["rho"] = env["density"]
    if "Nref" not in env and "N_axial" in env:
        env["Nref"] = env["N_axial"]
    if "k" not in env:
        if "kx" in env and "ky" in env:
            try:
                if float(env["kx"]) == float(env["ky"]):
                    env["k"] = env["kx"]
            except (TypeError, ValueError):
                pass
        elif "kx" in env:
            env["k"] = env["kx"]
        elif "ky" in env:
            env["k"] = env["ky"]
    if "L" not in env and "nodes" in env:
        nodes_arr = np.asarray(env["nodes"], dtype=float)
        if nodes_arr.ndim == 1:
            env["L"] = float(nodes_arr.max() - nodes_arr.min())
        elif nodes_arr.shape[1] >= 1:
            env["L"] = float(nodes_arr[:, 0].max() - nodes_arr[:, 0].min())
    if "nel" not in env and "nodes" in env:
        nodes_arr = np.asarray(env["nodes"])
        if nodes_arr.ndim >= 1 and nodes_arr.shape[0] >= 2:
            env["nel"] = int(nodes_arr.shape[0] - 1)

    _canonicalize_fix_entries(env, meta)
    _extract_element_properties(env, meta)
    _ensure_rect_mesh(env)

    code = meta.get("code", "")
    if code.startswith("THERMAL.") or "HEAT" in code.upper():
        _apply_temperature_dirichlet(env)

    return env


def _ensure_rect_mesh(env: Dict[str, Any]) -> None:
    if "nodes" in env or not all(k in env for k in ("Lx", "Ly", "nx", "ny")):
        return
    try:
        from .solvers.mesh_rect import rect_mesh
    except Exception:
        return
    Lx = float(env.get("Lx")); Ly = float(env.get("Ly"))
    nx = int(env.get("nx")); ny = int(env.get("ny"))
    nodes, elems = rect_mesh(Lx, Ly, nx, ny)
    env.setdefault("nodes", nodes)
    env.setdefault("elems", elems)
    env.setdefault("etype", "Q4")


def _canonicalize_fix_entries(env: Dict[str, Any], meta: Dict[str, Any]) -> None:
    entries = env.get("fix")
    if not entries:
        return
    categories = tuple(meta.get("categories") or ())
    normalized: List[List[Any]] = []
    for entry in entries:
        normalized.extend(_expand_fix_entry(entry, categories))
    if normalized:
        env["fix"] = normalized


def _expand_fix_entry(entry: Any, categories: Sequence[str]) -> List[List[Any]]:
    if entry is None:
        return []
    parts = list(entry)
    if not parts:
        return []
    node = int(parts[0])
    rest = parts[1:]
    if not rest:
        return []
    first = rest[0]
    if isinstance(first, str):
        if len(rest) == 1:
            return [[node, first, 0.0]]
        return [[node, first, rest[1]]]
    numeric = all(_is_number(v) for v in rest)
    if numeric and len(rest) >= 2:
        dofs = _infer_support_dofs(categories, len(rest))
        if dofs:
            out: List[List[Any]] = []
            for flag, dof in zip(rest, dofs):
                if float(flag):
                    out.append([node, dof, 0.0])
            return out
    if numeric and len(rest) == 1:
        return [[node, float(rest[0])]]
    return []


def _infer_support_dofs(categories: Sequence[str], count: int) -> Optional[Tuple[str, ...]]:
    cats = {str(cat).lower() for cat in categories}
    if count == 2:
        return ("ux", "uy")
    if count == 3:
        if "frame2d" in cats:
            return ("ux", "uy", "rz")
        if "truss3d" in cats or "solid" in cats:
            return ("ux", "uy", "uz")
    return None


def _is_number(value: Any) -> bool:
    try:
        float(value)
        return True
    except Exception:
        return False


def _extract_element_properties(env: Dict[str, Any], meta: Dict[str, Any]) -> None:
    categories = tuple(meta.get("categories") or ())
    allowed = {"truss2d", "frame2d", "beam", "bar1d", "contact", "truss3d"}
    if not any(cat in allowed for cat in categories):
        return
    elems = env.get("elems")
    if elems is None:
        return
    arr = np.asarray(elems)
    if arr.ndim != 2 or arr.shape[1] <= 3:
        return
    base = arr[:, :2].astype(int)
    extras = arr[:, 2:]
    env["elems"] = base.tolist()
    names = ["A", "E", "I", "prop3", "prop4"]
    for idx in range(min(len(names), extras.shape[1])):
        key = names[idx]
        if key in env:
            continue
        values = extras[:, idx]
        if np.allclose(values, values[0]):
            env[key] = float(values[0])
        else:
            env[key] = values.tolist()


def _apply_temperature_dirichlet(env: Dict[str, Any]) -> None:
    nodes = env.get("nodes")
    if nodes is None:
        return
    nodes = np.asarray(nodes, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] < 2:
        return
    tol = max(float(env.get("bc_tol", 1e-9)), 1e-9)
    fix_entries = [list(item) for item in env.get("fix", [])]
    existing_nodes = {int(item[0]) for item in fix_entries if len(item) >= 2}
    minx, maxx = nodes[:, 0].min(), nodes[:, 0].max()
    miny, maxy = nodes[:, 1].min(), nodes[:, 1].max()

    def add(mask, value):
        if value is None:
            return
        value = float(value)
        for idx in np.flatnonzero(mask):
            if idx in existing_nodes:
                continue
            fix_entries.append([int(idx), value])

    add(np.abs(nodes[:, 0] - minx) < tol, env.get("T_left"))
    add(np.abs(nodes[:, 0] - maxx) < tol, env.get("T_right"))
    add(np.abs(nodes[:, 1] - miny) < tol, env.get("T_bottom"))
    add(np.abs(nodes[:, 1] - maxy) < tol, env.get("T_top"))
    if fix_entries:
        env["fix"] = fix_entries


def _load_static_solver(code: str) -> Dict[str, Any] | None:
    spec = _STATIC_SOLVER_SPECS.get(code)
    if spec is None:
        return None
    module_name, func_name, required = spec
    module = import_module(f".solvers.{module_name}", package=__package__)
    func = getattr(module, func_name)
    entry = _entry(func, list(required))
    _SOLVER_REGISTRY[code] = entry
    return entry


def _load_hint_solver(
    code: str,
    hint: Optional[Tuple[str, Optional[str], Optional[List[str]]]],
    categories: Sequence[str],
    objective: str,
) -> Optional[Dict[str, Any]]:
    if not hint:
        return None
    module_name, func_name, required = hint
    module = import_module(f".solvers.{module_name}", package=__package__)
    func = None
    if func_name:
        func = getattr(module, func_name, None)
    if func is None:
        func = _resolve_solver_function(module, module_name, categories, objective)
    if func is None:
        fallback = getattr(module, f"solve_{module_name}", None)
        if callable(fallback):
            func = fallback
    if func is None:
        candidates = [
            getattr(module, name)
            for name in dir(module)
            if name.startswith("solve_") and callable(getattr(module, name))
        ]
        if candidates:
            func = candidates[0]
    if func is None:
        return None
    inputs = required
    if inputs is None:
        inputs = getattr(module, "REQUIRED_INPUTS", None)
        if inputs is None:
            inputs = getattr(func, "required_inputs", None)
    entry = {"func": func, "required": list(inputs or [])}
    _SOLVER_REGISTRY[code] = entry
    return entry


def _to_snake(token: str) -> str:
    text = token or ""
    text = text.replace("-", "_").replace(" ", "_").replace(".", "_")
    text = re.sub(r"([a-zA-Z])(\d)", r"\1_\2", text)
    text = re.sub(r"(\d)([a-zA-Z])", r"\1_\2", text)
    text = re.sub(r"(?<!^)(?=[A-Z])", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_").lower()


def _auto_discover_solver(
    categories: Sequence[str],
    objective: str,
    raw_categories: Sequence[str],
    raw_objective: str,
) -> Dict[str, Any] | None:
    tokens = [tok for tok in categories if tok]
    objective = (objective or "").strip()
    raw_cat_snake = [_to_snake(tok) for tok in raw_categories if tok]
    objective_snake = _to_snake(raw_objective)
    module_candidates: List[str] = []

    for obj in filter(None, [objective_snake, objective]):
        if raw_cat_snake:
            module_candidates.append("_".join(raw_cat_snake + [obj]))
        if tokens:
            module_candidates.append("_".join(tokens + [obj]))
    if raw_cat_snake:
        module_candidates.append("_".join(raw_cat_snake))
    if tokens:
        module_candidates.append("_".join(tokens))
    module_candidates.extend(raw_cat_snake)
    module_candidates.extend(tokens)
    module_candidates.extend(filter(None, [objective_snake, objective]))

    tried = set()
    for module_name in module_candidates:
        if not module_name or module_name in tried:
            continue
        tried.add(module_name)
        try:
            module = import_module(f".solvers.{module_name}", package=__package__)
        except ModuleNotFoundError:
            continue
        func = _resolve_solver_function(module, module_name, tokens, objective)
        if func is None:
            continue
        required = getattr(module, "REQUIRED_INPUTS", None)
        if required is None:
            required = getattr(func, "required_inputs", None)
        req_list = list(required or [])
        return {"func": func, "required": req_list}
    return None


def _resolve_solver_function(module, module_name: str, tokens: Sequence[str], objective: str):
    base = "_".join(tokens)
    candidates = []
    if module_name:
        candidates.append(f"solve_{module_name}")
    if base and objective:
        candidates.append(f"solve_{base}_{objective}")
    if base:
        candidates.append(f"solve_{base}")
    if objective:
        candidates.append(f"solve_{objective}")
    for name in dict.fromkeys(candidates):
        func = getattr(module, name, None)
        if callable(func):
            return func
    return None


def _get_solver_entry(
    code: str,
    categories: Sequence[str],
    objective: str,
    raw_categories: Sequence[str],
    raw_objective: str,
) -> Dict[str, Any]:
    entry = _SOLVER_REGISTRY.get(code)
    if entry is not None:
        return entry
    loaded = _load_static_solver(code)
    if loaded is not None:
        return loaded
    hint_entry = _load_hint_solver(
        code, _get_solver_hints().get(code), categories, objective
    )
    if hint_entry is not None:
        return hint_entry
    discovered = _auto_discover_solver(
        categories,
        objective,
        raw_categories,
        raw_objective,
    )
    if discovered is None:
        raise NotImplementedError(
            f'Solver not implemented for PROBLEM code "{code}".'
        )
    _SOLVER_REGISTRY[code] = discovered
    return discovered


def dispatch(problem: Dict[str, Any], env: Dict[str, Any]) -> Dict[str, Any]:
    meta = _problem_metadata(problem or {})
    code = meta["code"]
    categories = meta.get("categories") or ()
    objective = meta.get("objective") or ""
    raw_categories = meta.get("raw_categories") or ()
    raw_objective = meta.get("raw_objective") or ""
    env = _normalize_environment(env, meta)
    entry = _get_solver_entry(code, categories, objective, raw_categories, raw_objective)
    required = entry["required"]
    if required:
        missing = [name for name in required if name not in env]
        if missing:
            raise KeyError(f"Inputs {missing} required for {code}")
    return entry["func"](env)
