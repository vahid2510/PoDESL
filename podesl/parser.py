from __future__ import annotations

import codeop
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from .problem_codes import parse_problem_title, problem_code
from .textutils import split_semicolons, strip_inline_comment

_HEADERS = ("GIVEN", "EQUATIONS", "SOLVE", "REPORT")


class DSLParseError(Exception):
    """Raised when DSL source cannot be parsed."""


@dataclass
class _Statement:
    kind: str  # "assign" or "command"
    name: str
    expr: Optional[str] = None
    args: Optional[List[str]] = None
    line: int = 0


class _Parser:
    def __init__(self, source: str, filename: str) -> None:
        self.source = source
        self.filename = filename
        self.lines = source.splitlines()
        self.total = len(self.lines)
        self.index = 0  # 0-based
        self.problem: Optional[Dict[str, Any]] = None
        self.blocks: List[Dict[str, Any]] = []

    def parse(self) -> Dict[str, Any]:
        self._parse_problem()
        while self.index < self.total:
            line = self._current_line()
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                self.index += 1
                continue
            upper = stripped.upper()
            if upper in _HEADERS:
                self._parse_block(upper)
                continue
            raise DSLParseError(
                f"{self.filename}:{self.index + 1}: Unexpected statement outside block: {line.strip()}"
            )
        if not self.problem:
            raise DSLParseError(f"{self.filename}: Missing PROBLEM definition.")
        return {"problem": self.problem, "blocks": self.blocks}

    def _parse_problem(self) -> None:
        while self.index < self.total:
            line = self._current_line()
            stripped = line.strip()
            if not stripped:
                self.index += 1
                continue
            if stripped.startswith("#"):
                self.index += 1
                continue
            if stripped.upper().startswith("PROBLEM"):
                title = self._extract_title(stripped)
                descriptor = parse_problem_title(title)
                canonical_code = problem_code(title)
                self.problem = {
                    "type": "problem",
                    "title": title,
                    "domain": descriptor.domain,
                    "ptype": descriptor.primary_category,
                    "categories": list(descriptor.categories),
                    "raw_categories": list(descriptor.raw_categories),
                    "objective": descriptor.objective,
                    "raw_objective": descriptor.raw_objective,
                    "code": canonical_code,
                }
                self.index += 1
                return
            raise DSLParseError(
                f"{self.filename}:{self.index + 1}: Expected PROBLEM definition."
            )
        raise DSLParseError(f"{self.filename}: Empty source.")

    def _extract_title(self, line: str) -> str:
        rest = line[len("PROBLEM") :].strip()
        if not rest:
            raise DSLParseError(
                f"{self.filename}:{self.index + 1}: PROBLEM title is missing."
            )
        if rest[0] in "\"'":
            quote = rest[0]
            if rest.count(quote) < 2 or not rest.endswith(quote):
                raise DSLParseError(
                    f"{self.filename}:{self.index + 1}: Unterminated PROBLEM title."
                )
            return rest[1:-1]
        return rest

    def _parse_block(self, header: str) -> None:
        block = {"type": header.lower(), "stmts": []}
        self.blocks.append(block)
        self.index += 1
        while self.index < self.total:
            line = self._current_line()
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                self.index += 1
                continue
            upper = stripped.upper()
            if upper in _HEADERS and line.lstrip() == line:
                # Next block begins (no indentation)
                break
            if line.lstrip() == line:
                raise DSLParseError(
                    f"{self.filename}:{self.index + 1}: Statements inside {header} must be indented."
                )
            self._parse_statements_into_block(block, line)
            self.index += 1

    def _parse_statements_into_block(self, block: Dict[str, Any], line: str) -> None:
        fragments = split_semicolons(line)
        for fragment in fragments:
            cleaned = strip_inline_comment(fragment)
            for piece in _split_inline_assignments(cleaned):
                text = piece.strip()
                if not text:
                    continue
                stmt = self._parse_statement(text)
                if stmt:
                    block["stmts"].append(
                        {
                            "kind": stmt.kind,
                            "name": stmt.name,
                            "value": stmt.expr,
                            "args": stmt.args or [],
                            "line": self.index + 1,
                        }
                    )

    def _parse_statement(self, text: str) -> Optional[_Statement]:
        if "=" in text:
            name, expr = self._split_assignment(text)
            if name:
                value = self._collect_expression(expr)
                return _Statement(kind="assign", name=name, expr=value)
        parts = text.split(None, 1)
        if not parts:
            return None
        name = parts[0]
        arg_str = parts[1] if len(parts) > 1 else ""
        args = self._split_args(arg_str)
        return _Statement(kind="command", name=name, args=args)

    def _split_assignment(self, text: str) -> tuple[Optional[str], str]:
        name = []
        idx = 0
        while idx < len(text) and not text[idx].isspace() and text[idx] != "=":
            name.append(text[idx])
            idx += 1
        if not name:
            return None, text
        remaining = text[idx:].lstrip()
        if not remaining.startswith("="):
            return None, text
        expr = remaining[1:].strip()
        return "".join(name), expr

    def _collect_expression(self, initial: str) -> str:
        expr_lines: List[str] = []
        current_line = self.index
        if initial:
            expr_lines.append(initial)
        while True:
            joined = "\n".join(expr_lines).strip()
            if joined:
                try:
                    code = codeop.compile_command(
                        joined, filename=self.filename, symbol="eval"
                    )
                except SyntaxError as exc:  # pragma: no cover - precise feedback
                    raise DSLParseError(
                        f"{self.filename}:{self.index + 1}: Invalid expression: {exc.msg}"
                    ) from exc
                if code is not None:
                    return joined
            current_line += 1
            if current_line >= self.total:
                raise DSLParseError(
                    f"{self.filename}:{self.index + 1}: Unterminated expression."
                )
            next_line = strip_inline_comment(self.lines[current_line]).strip()
            expr_lines.append(next_line)
            self.index = current_line

    def _split_args(self, text: str) -> List[str]:
        text = strip_inline_comment(text or "").strip()
        if not text:
            return []
        args: List[str] = []
        buf: List[str] = []
        depth = 0
        quote: Optional[str] = None
        escape = False
        for ch in text:
            if escape:
                buf.append(ch)
                escape = False
                continue
            if ch == "\\" and quote:
                buf.append(ch)
                escape = True
                continue
            if ch in ("'", '"'):
                if quote == ch:
                    quote = None
                elif quote is None:
                    quote = ch
                buf.append(ch)
                continue
            if quote is None:
                if ch in "([{":
                    depth += 1
                elif ch in ")]}":
                    depth = max(0, depth - 1)
                elif ch.isspace() and depth == 0:
                    if buf:
                        args.append("".join(buf).strip())
                        buf = []
                    continue
            buf.append(ch)
        if buf:
            args.append("".join(buf).strip())
        return args

    def _current_line(self) -> str:
        return self.lines[self.index]


def parse(src: str, filename: str = "<stdin>") -> Dict[str, Any]:
    code = src if isinstance(src, str) else src.read()
    parser = _Parser(code, filename)
    return parser.parse()


def _split_inline_assignments(text: str) -> List[str]:
    if not text:
        return []
    parts: List[str] = []
    buf: List[str] = []
    quote: Optional[str] = None
    depth = 0
    idx = 0
    length = len(text)

    def find_assignment(start: int) -> Optional[int]:
        pos = start
        while pos < length and text[pos].isspace():
            pos += 1
        if pos >= length:
            return None
        if not (text[pos].isalpha() or text[pos] == "_"):
            return None
        end = pos + 1
        while end < length and (text[end].isalnum() or text[end] == "_"):
            end += 1
        if end < length and text[end] == "=":
            return pos
        return None

    while idx < length:
        ch = text[idx]
        if quote:
            if ch == "\\":
                buf.append(ch)
                idx += 1
                if idx < length:
                    buf.append(text[idx])
                    idx += 1
                continue
            if ch == quote:
                quote = None
            buf.append(ch)
            idx += 1
            continue
        if ch in ("'", '"'):
            quote = ch
            buf.append(ch)
            idx += 1
            continue
        if ch in "([{":
            depth += 1
            buf.append(ch)
            idx += 1
            continue
        if ch in ")]}":
            depth = max(0, depth - 1)
            buf.append(ch)
            idx += 1
            continue
        if depth == 0 and ch.isspace():
            next_pos = find_assignment(idx + 1)
            if next_pos is not None:
                segment = "".join(buf).strip()
                if segment:
                    parts.append(segment)
                buf = []
                idx = next_pos
                continue
        buf.append(ch)
        idx += 1

    if buf:
        segment = "".join(buf).strip()
        if segment:
            parts.append(segment)

    return parts
