from __future__ import annotations
from typing import List
from .helpdb import HELP_INDEX
from .problem_codes import is_supported, problem_code

def list_problems() -> List[str]:
    return [meta["title"] for _, meta in HELP_INDEX.items()]

def help_problem(title_or_code: str) -> str:
    code = title_or_code
    if code not in HELP_INDEX:
        if is_supported(title_or_code):
            code = problem_code(title_or_code)
    if code not in HELP_INDEX:
        raise KeyError(f"No help available for: {title_or_code}")
    meta = HELP_INDEX[code]
    lines = []
    lines.append(f"PROBLEM: {meta['title']}  [{code}]")
    lines.append(f"Summary: {meta['summary']}")
    lines.append("")
    lines.append("Required inputs:")
    for k in meta["inputs_required"]: lines.append(f"  - {k}")
    if meta.get("inputs_optional"):
        lines.append("Optional inputs:")
        for k in meta["inputs_optional"]: lines.append(f"  - {k}")
    lines.append("Outputs:")
    for k in meta["outputs"]: lines.append(f"  - {k}")
    lines.append("")
    lines.append("Example:")
    lines.append(meta["example"])
    return "\n".join(lines)

def help_all() -> str:
    out = []
    for code in sorted(HELP_INDEX.keys()):
        out.append(help_problem(code))
        out.append("-"*60)
    return "\n".join(out)
