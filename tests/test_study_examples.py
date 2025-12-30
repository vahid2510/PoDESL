from __future__ import annotations

from pathlib import Path

from podesl.cli import run_file

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "examples"


def _run_with_outputs(name: str, outputs: list[str]) -> None:
    # Clean up any stale artifacts
    for rel in outputs:
        target = ROOT / rel
        if target.exists():
            target.unlink()

    run_file(EXAMPLES_DIR / name, emit_python=False, save_python=None, dry_run=False)

    for rel in outputs:
        target = ROOT / rel
        assert target.exists(), f"Expected output {rel} from {name}"
        target.unlink()


def test_example_study_montecarlo_runs():
    _run_with_outputs(
        "study_montecarlo_bar.dsl",
        ["study_montecarlo_mean.csv", "study_montecarlo_std.csv"],
    )


def test_example_study_optimize_runs():
    _run_with_outputs(
        "study_optimize_bar.dsl",
        ["study_optimize_x_opt.csv", "study_optimize_history.csv"],
    )


def test_example_study_scenarios_runs():
    _run_with_outputs(
        "study_scenarios_bar.dsl",
        ["study_scenarios_tip_disp.csv"],
    )


def test_example_study_fusion_runs():
    _run_with_outputs("study_fusion_steady.dsl", [])
    _run_with_outputs(
        "study_fusion_transient.dsl",
        ["study_fusion_transient_times.csv", "study_fusion_transient_tip.csv"],
    )
