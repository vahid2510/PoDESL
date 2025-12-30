from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional
import runpy

# Ensure package imports work when executed as a standalone script or frozen binary.
if __package__ in (None, ""):
    _ROOT = Path(__file__).resolve().parents[1]
    if str(_ROOT) not in sys.path:
        sys.path.insert(0, str(_ROOT))

from podesl.parser import parse
from podesl.transpile import transpile
from podesl.help import list_problems, help_problem, help_all

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def _emit_line(text: str = "") -> None:
    try:
        print(text)
    except UnicodeEncodeError:
        encoding = sys.stdout.encoding or "utf-8"
        sys.stdout.buffer.write(text.encode(encoding, errors="replace"))
        sys.stdout.buffer.write(b"\n")


def _load_source(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8").lstrip("\ufeff")
    except FileNotFoundError:
        raise SystemExit(f'Input file not found: "{path}"')


def _execute_transpiled(code: str, origin: str) -> None:
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", suffix=".py", delete=False
        ) as tmp:
            tmp.write(code)
            tmp_path = Path(tmp.name)
        env = os.environ.copy()
        existing = env.get("PYTHONPATH")
        path_entries = [str(PACKAGE_ROOT)]
        if existing:
            path_entries.append(existing)
        env["PYTHONPATH"] = os.pathsep.join(path_entries)
        if getattr(sys, "frozen", False):
            saved_argv = sys.argv[:]
            sys.argv = [str(tmp_path)]
            saved_py_path = os.environ.get("PYTHONPATH")
            os.environ["PYTHONPATH"] = env["PYTHONPATH"]
            try:
                runpy.run_path(str(tmp_path), run_name="__main__")
            finally:
                if saved_py_path is None:
                    os.environ.pop("PYTHONPATH", None)
                else:
                    os.environ["PYTHONPATH"] = saved_py_path
                sys.argv = saved_argv
        else:
            subprocess.run([sys.executable, str(tmp_path)], check=True, env=env)
    except subprocess.CalledProcessError as exc:
        msg = (
            f"Execution failed for {origin} (exit code {exc.returncode}). "
            "Re-run with --emit-python to inspect the generated script."
        )
        raise SystemExit(msg) from exc
    finally:
        if tmp_path:
            try:
                tmp_path.unlink()
            except FileNotFoundError:
                pass


def run_file(path: Path, emit_python: bool, save_python: Optional[Path], dry_run: bool) -> None:
    src = _load_source(path)
    try:
        ast = parse(src, filename=str(path))
    except Exception as exc:
        raise SystemExit(f"Failed to parse {path}: {exc}") from exc

    try:
        code = transpile(ast)
    except Exception as exc:
        raise SystemExit(f"Failed to transpile {path}: {exc}") from exc

    if emit_python:
        print(code)

    if save_python:
        save_python.write_text(code, encoding="utf-8")

    if not dry_run:
        _execute_transpiled(code, str(path))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="podesl command line interface")
    parser.add_argument("input", nargs="?", help="Path to a .dsl problem file")
    parser.add_argument("--list-problems", action="store_true", dest="list_problems",
                        help="List supported PROBLEM definitions")
    parser.add_argument("--help-problem", metavar="TITLE_OR_CODE",
                        help="Show help for a given PROBLEM title or code")
    parser.add_argument("--help-all", action="store_true",
                        help="Dump help for all supported problems")
    parser.add_argument("--emit-python", action="store_true",
                        help="Print transpiled Python to stdout")
    parser.add_argument("--save-python", metavar="PATH",
                        help="Write transpiled Python to the provided path")
    parser.add_argument("--dry-run", action="store_true",
                        help="Skip execution after transpiling")
    parser.add_argument("--export-abaqus", metavar="PATH",
                        help="Write an Abaqus .inp deck for the provided DSL file")
    parser.add_argument("--abaqus-job", metavar="NAME",
                        help="Override the Abaqus job name when exporting")
    parser.add_argument("--ide", action="store_true",
                        help="Launch the interactive web IDE instead of executing a single file")
    parser.add_argument("--ide-host", default="127.0.0.1",
                        help="Host/IP for the IDE server (default: 127.0.0.1)")
    parser.add_argument("--ide-port", type=int, default=8765,
                        help="Port for the IDE server (default: 8765)")
    parser.add_argument("--ide-no-browser", action="store_true",
                        help="Do not automatically open a browser window for the IDE")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if getattr(args, "ide", False):
        from podesl.ide.server import run_ide_server

        run_ide_server(
            host=args.ide_host,
            port=args.ide_port,
            open_browser=not args.ide_no_browser,
        )
        return 0

    if args.list_problems:
        for title in sorted(list_problems()):
            _emit_line(title)
        return 0

    if args.help_all:
        for line in help_all().splitlines():
            _emit_line(line)
        return 0

    if args.help_problem:
        for line in help_problem(args.help_problem).splitlines():
            _emit_line(line)
        return 0

    if args.export_abaqus and not args.input:
        parser.error("--export-abaqus requires an input DSL file.")

    if not args.input:
        parser.print_help()
        return 1

    save_path = Path(args.save_python).resolve() if args.save_python else None
    if args.export_abaqus:
        from podesl.abaqus_export import export_abaqus_input

        source_text = Path(args.input).read_text(encoding="utf-8")
        metadata = {}
        if args.abaqus_job:
            metadata["job_name"] = args.abaqus_job
        abaqus_path = Path(args.export_abaqus).resolve()
        info = export_abaqus_input(source_text, args.input, abaqus_path, metadata=metadata)
        _emit_line(
            f"Abaqus input saved to {info.path} "
            f"(nodes={info.node_count}, elements={info.element_count})"
        )

    run_file(Path(args.input), args.emit_python, save_path, args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
