from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import List, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from podesl.cli import run_file


def _example_files(root: Path) -> List[Path]:
    examples_dir = root / "examples"
    return sorted(examples_dir.glob("*.dsl"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parse/transpile every DSL example to catch regressions."
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Run each transpiled script instead of a dry-run (can be slow)",
    )
    args = parser.parse_args()

    examples = _example_files(ROOT)
    failures: List[Tuple[Path, str]] = []
    for dsl_path in examples:
        try:
            run_file(dsl_path, emit_python=False, save_python=None, dry_run=not args.execute)
        except SystemExit as exc:
            message = str(exc) or f"exit code {exc.code}"
            failures.append((dsl_path, message))
        except Exception as exc:  # pragma: no cover - defensive
            failures.append((dsl_path, f"unexpected error: {exc}"))

    if failures:
        print("Smoke test failures:\n")
        for path, message in failures:
            print(f"- {path}: {message}")
        return 1

    print(f"Validated {len(examples)} example(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
