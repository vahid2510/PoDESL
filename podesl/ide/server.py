from __future__ import annotations

import json
import os
import re
import threading
import time
import urllib.parse
import webbrowser
from collections import deque
from functools import partial
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
from uuid import uuid4

import numpy as np

from ..abaqus_export import export_abaqus_input
from ..execution import execute_source, evaluate_last_expression, get_last_run
from ..help import help_problem, list_problems
from ..parser import DSLParseError

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
STATIC_ROOT = Path(__file__).resolve().with_name("static")
EXAMPLES_ROOT = PACKAGE_ROOT / "examples"
DEFAULT_BUILD_DIR = PACKAGE_ROOT / "build" / "ide_outputs"

DSL_HEADERS = ("PROBLEM", "GIVEN", "SOLVE", "EQUATIONS", "REPORT")
EXAMPLE_PATTERNS = ("*.dsl", "*.md", "*.rst", "*.txt")
KEYWORDS = [
    "PROBLEM",
    "GIVEN",
    "SOLVE",
    "EQUATIONS",
    "REPORT",
    "REPORT print",
    "REPORT export",
    "REPORT plot",
    "help",
    "linspace",
    "zeros",
    "ones",
    "array",
    "diag",
]

RUN_HISTORY: deque[Dict[str, Any]] = deque(maxlen=32)


def _json_serialize(value: Any, depth: int = 0) -> Any:
    if depth > 4:
        return "... "
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _json_serialize(v, depth + 1) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_serialize(item, depth + 1) for item in value]
    if hasattr(value, "tolist"):
        try:
            return np.asarray(value).tolist()
        except Exception:
            pass
    return repr(value)


def _format_dsl(source: str, indent: int = 2) -> str:
    lines = source.splitlines()
    formatted: List[str] = []
    level = 0
    for raw in lines:
        stripped = raw.strip()
        if not stripped:
            formatted.append("")
            continue
        upper = stripped.upper()
        if upper.startswith("PROBLEM"):
            formatted.append(stripped)
            level = 0
            continue
        if upper in DSL_HEADERS:
            formatted.append(upper.title())
            level = 1
            continue
        indent_prefix = " " * (indent * level)
        formatted.append(f"{indent_prefix}{stripped}")
    return "\n".join(formatted)


def _list_examples() -> List[Dict[str, Any]]:
    entries: List[Dict[str, Any]] = []
    if not EXAMPLES_ROOT.exists():
        return entries
    seen: Dict[str, Dict[str, Any]] = {}
    for pattern in EXAMPLE_PATTERNS:
        for path in EXAMPLES_ROOT.rglob(pattern):
            try:
                rel = path.relative_to(PACKAGE_ROOT)
            except ValueError:
                rel = path
            key = str(rel).replace(os.sep, "/")
            if key in seen:
                continue
            stat = path.stat()
            suffix = path.suffix.lower()
            kind = "doc" if suffix in {".md", ".rst", ".txt"} else "dsl"
            seen[key] = {
                "path": key,
                "name": path.stem,
                "size": stat.st_size,
                "modified": stat.st_mtime,
                "kind": kind,
            }
    entries.extend(seen.values())
    entries.sort(key=lambda item: item["path"].lower())
    return entries


def _ensure_series(value: Any) -> List[float]:
    arr = np.asarray(value, dtype=float).ravel()
    return arr.tolist()


def _record_history(entry: Dict[str, Any]) -> Dict[str, Any]:
    payload = dict(entry)
    payload.setdefault("id", uuid4().hex)
    payload.setdefault("timestamp", time.time())
    RUN_HISTORY.appendleft(payload)
    return payload


def _diagnostics_from_exception(exc: Exception) -> List[Dict[str, Any]]:
    diags: List[Dict[str, Any]] = []
    text = str(exc)
    match = re.search(r":(\d+):", text)
    if match:
        line = int(match.group(1))
        diags.append({"line": line, "column": 0, "message": text})
        return diags
    match = re.search(r"line (\d+)", text, re.IGNORECASE)
    if match:
        line = int(match.group(1))
        diags.append({"line": line, "column": 0, "message": text})
    return diags


def _resolve_workspace_path(rel: str) -> Path:
    rel = rel.lstrip("/\\")
    candidate = (PACKAGE_ROOT / rel).resolve()
    if PACKAGE_ROOT not in candidate.parents and candidate != PACKAGE_ROOT:
        raise PermissionError("Path escapes workspace.")
    return candidate


def _gather_identifiers(source: str) -> List[str]:
    identifiers: List[str] = []
    pattern = re.compile(r"^\s*([A-Za-z_]\w*)\s*=", re.MULTILINE)
    for match in pattern.finditer(source):
        name = match.group(1)
        if name not in identifiers:
            identifiers.append(name)
    return identifiers


def _current_token(source: str, cursor: int) -> str:
    cursor = max(0, min(len(source), cursor))
    start = cursor
    while start > 0 and (source[start - 1].isalnum() or source[start - 1] in ("_", ":", ".")):
        start -= 1
    end = cursor
    while end < len(source) and (source[end].isalnum() or source[end] in ("_", ":", ".")):
        end += 1
    return source[start:end]


def _autocomplete(
    source: str,
    cursor: int,
    problems: Iterable[str],
) -> Dict[str, Any]:
    token = _current_token(source, cursor)
    token_lower = token.lower()
    suggestions: List[Dict[str, Any]] = []

    def append(label: str, insert: Optional[str] = None, kind: str = "keyword") -> None:
        if token_lower and not label.lower().startswith(token_lower):
            return
        suggestions.append(
            {
                "label": label,
                "insertText": insert or label,
                "kind": kind,
            }
        )

    for kw in KEYWORDS:
        append(kw)

    for pid in problems:
        append(pid, insert=f'"{pid}"', kind="problem")

    for ident in _gather_identifiers(source):
        append(ident, kind="variable")

    return {"token": token, "suggestions": suggestions[:50]}


class _IDEResponseMixin:
    def _send_json(self, payload: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status.value)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_json_body(self) -> Dict[str, Any]:
        length = int(self.headers.get("Content-Length", 0))
        data = self.rfile.read(length) if length else b""
        if not data:
            return {}
        try:
            return json.loads(data.decode("utf-8"))
        except json.JSONDecodeError:
            return {}

    def _json_error(self, message: str, status: HTTPStatus = HTTPStatus.BAD_REQUEST) -> None:
        self._send_json({"error": message}, status=status)


class IDEServerHandler(_IDEResponseMixin, SimpleHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def __init__(self, *args, directory: Optional[str] = None, **kwargs):
        if directory is None:
            directory = str(STATIC_ROOT)
        super().__init__(*args, directory=directory, **kwargs)

    def log_message(self, format: str, *args: Any) -> None:
        # Quieter logging to avoid spamming the console during editing.
        message = "%s - - [%s] %s\n" % (
            self.address_string(),
            self.log_date_time_string(),
            format % args,
        )
        print(message.rstrip())

    def do_GET(self) -> None:  # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path.startswith("/api/"):
            self._handle_api_get(parsed)
            return
        if parsed.path == "/":
            self.path = "/index.html"
        return super().do_GET()

    def do_POST(self) -> None:  # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        if not parsed.path.startswith("/api/"):
            return self.send_error(HTTPStatus.NOT_FOUND.value)
        self._handle_api_post(parsed)

    def _handle_api_get(self, parsed: urllib.parse.ParseResult) -> None:
        route = parsed.path
        query = urllib.parse.parse_qs(parsed.query)
        if route == "/api/examples":
            self._send_json({"examples": _list_examples()})
            return
        if route == "/api/file":
            rel = query.get("path", [None])[0]
            if not rel:
                return self._json_error("path query parameter missing")
            try:
                target = _resolve_workspace_path(rel)
            except PermissionError as exc:
                return self._json_error(str(exc), status=HTTPStatus.FORBIDDEN)
            if not target.exists():
                return self._json_error("File not found.", status=HTTPStatus.NOT_FOUND)
            text = target.read_text(encoding="utf-8")
            self._send_json({"path": rel, "content": text})
            return
        if route == "/api/help":
            term = query.get("query", [""])[0].strip()
            problems = list_problems()
            if not term:
                matches = problems
                detail = None
            else:
                try:
                    detail = help_problem(term)
                    matches = [term]
                except Exception:
                    detail = None
                    matches = [p for p in problems if term.lower() in p.lower()]
            self._send_json({"query": term, "matches": matches[:40], "detail": detail})
            return
        if route == "/api/problems":
            self._send_json({"problems": list_problems()})
            return
        if route == "/api/history":
            self._send_json({"history": list(RUN_HISTORY)})
            return
        if route == "/api/status":
            self._send_json(
                {
                    "workspace": str(PACKAGE_ROOT),
                    "examples": len(_list_examples()),
                    "time": time.time(),
                }
            )
            return
        self._send_json({"error": "Unknown endpoint."}, status=HTTPStatus.NOT_FOUND)

    def _handle_api_post(self, parsed: urllib.parse.ParseResult) -> None:
        route = parsed.path
        payload = self._read_json_body()
        if route == "/api/run":
            source = payload.get("source", "")
            filename = payload.get("filename") or "<browser>"
            try:
                result = execute_source(source, filename=filename)
            except Exception as exc:
                diagnostics = _diagnostics_from_exception(exc)
                _record_history(
                    {
                        "status": "error",
                        "path": filename,
                        "message": str(exc),
                        "diagnostics": diagnostics,
                    }
                )
                self._send_json(
                    {"error": str(exc), "diagnostics": diagnostics},
                    status=HTTPStatus.BAD_REQUEST,
                )
                return
            data = {
                "stdout": result.stdout,
                "stderr": result.stderr,
                "result": _json_serialize(result.result),
                "problem": result.problem,
                "env": _json_serialize(result.env),
            }
            entry = _record_history(
                {
                    "status": "ok",
                    "path": filename,
                    "message": "OK",
                    "problem": result.problem,
                }
            )
            data["history_id"] = entry["id"]
            self._send_json(data)
            return
        if route == "/api/format":
            source = payload.get("source", "")
            formatted = _format_dsl(source)
            self._send_json({"text": formatted})
            return
        if route == "/api/save":
            rel = payload.get("path")
            if not rel:
                return self._json_error("Missing path.")
            content = payload.get("source", "")
            try:
                target = _resolve_workspace_path(rel)
            except PermissionError as exc:
                return self._json_error(str(exc), status=HTTPStatus.FORBIDDEN)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content, encoding="utf-8")
            self._send_json({"path": rel, "bytes": len(content.encode("utf-8"))})
            return
        if route == "/api/autocomplete":
            source = payload.get("source", "")
            cursor = int(payload.get("cursor", len(source)))
            self._send_json(_autocomplete(source, cursor, list_problems()))
            return
        if route == "/api/export":
            source = payload.get("source", "")
            filename = payload.get("filename") or "<browser>"
            desired = payload.get("output_path") or ""
            metadata = payload.get("metadata") or {}
            output = (
                _resolve_workspace_path(desired)
                if desired
                else DEFAULT_BUILD_DIR / f"{metadata.get('job_name') or 'podesl_model'}.inp"
            )
            output.parent.mkdir(parents=True, exist_ok=True)
            try:
                info = export_abaqus_input(source, filename, output, metadata=metadata)
            except Exception as exc:
                return self._json_error(f"Abaqus export failed: {exc}")
            self._send_json(
                {
                    "path": str(info.path),
                    "job_name": info.job_name,
                    "node_count": info.node_count,
                    "element_count": info.element_count,
                }
            )
            return
        if route == "/api/chart":
            x_expr = payload.get("x_expr")
            y_expr = payload.get("y_expr")
            if not x_expr or not y_expr:
                return self._json_error("chart requires x_expr and y_expr")
            try:
                x_vals = evaluate_last_expression(x_expr)
                y_vals = evaluate_last_expression(y_expr)
                x_series = _ensure_series(x_vals)
                y_series = _ensure_series(y_vals)
            except Exception as exc:
                return self._json_error(f"Chart evaluation failed: {exc}")
            if len(x_series) != len(y_series):
                return self._json_error("Series length mismatch between X and Y.")
            self._send_json({"x": x_series, "y": y_series})
            return
        if route == "/api/report":
            fields = payload.get("fields") or []
            notes = payload.get("notes") or ""
            try:
                ctx = get_last_run()
            except Exception as exc:
                return self._json_error(str(exc))
            metrics: List[Dict[str, Any]] = []
            for field in fields:
                expr = field.get("expr")
                label = field.get("label") or expr
                if not expr:
                    continue
                try:
                    value = evaluate_last_expression(expr)
                except Exception as exc:
                    value = f"error: {exc}"
                metrics.append(
                    {
                        "label": label,
                        "expr": expr,
                        "value": _json_serialize(value),
                    }
                )
            summary = {
                "problem": ctx.get("problem"),
                "notes": notes,
                "metrics": metrics,
            }
            self._send_json(summary)
            return
        if route == "/api/history":
            if payload.get("action") == "clear":
                RUN_HISTORY.clear()
                self._send_json({"history": []})
            else:
                self._send_json({"history": list(RUN_HISTORY)})
            return
        self._send_json({"error": "Unknown endpoint."}, status=HTTPStatus.NOT_FOUND)


def run_ide_server(host: str = "127.0.0.1", port: int = 8765, open_browser: bool = True) -> None:
    handler_cls = partial(IDEServerHandler, directory=str(STATIC_ROOT))
    httpd = ThreadingHTTPServer((host, port), handler_cls)
    url = f"http://{host}:{port}/"
    print(f"PoDESL IDE running at {url}")

    if open_browser:
        threading.Thread(target=lambda: webbrowser.open(url), daemon=True).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("Shutting down IDE server...")
    finally:
        httpd.server_close()
