from __future__ import annotations

import json
import pathlib
import time
import traceback

import sys

ROOT_DIR = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from podesl.execution import execute_source

ROOT = ROOT_DIR / "examples"

def main() -> None:
    results = []
    start = time.time()
    count = 0
    for path in sorted(ROOT.rglob("*.dsl")):
        count += 1
        src = path.read_text(encoding="utf-8")
        entry = {"path": str(path).replace("\\", "/")}
        try:
            execute_source(src, filename=str(path))
            entry["status"] = "ok"
        except Exception as exc:  # noqa: BLE001
            entry["status"] = "error"
            entry["error"] = repr(exc)
            entry["traceback"] = traceback.format_exc()
        results.append(entry)
    duration = time.time() - start
    payload = {"count": count, "duration": duration, "results": results}
    pathlib.Path("example_audit.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    failures = sum(1 for r in results if r["status"] != "ok")
    print(f"Audited {count} DSL files in {duration:.1f}s; failures: {failures}")


if __name__ == "__main__":
    main()
