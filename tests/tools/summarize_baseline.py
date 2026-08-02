#!/usr/bin/env python3
"""Summarizes tests/baseline/ into a single report.

Part of issue #11 task 4 (ROADMAP.md Phase 0 baseline capture). Walks
tests/baseline/<case>/<solver>_<mode>_<flow>/ directories produced by
run_baseline.sh / run_baseline_batch.sh and writes a machine-readable
JSON summary plus a human-readable table.

Usage: summarize_baseline.py [baseline_dir] [--json PATH]
"""

import argparse
import json
import sys
from pathlib import Path


def collect(baseline_dir: Path) -> list[dict]:
    rows = []
    for case_dir in sorted(baseline_dir.iterdir()):
        if not case_dir.is_dir():
            continue
        for run_dir in sorted(case_dir.iterdir()):
            verdict_file = run_dir / "verdict.txt"
            wall_file = run_dir / "wallclock_seconds.txt"
            exit_file = run_dir / "exit_status.txt"
            if not verdict_file.exists():
                continue
            solver_mode_flow = run_dir.name
            rows.append(
                {
                    "case": case_dir.name,
                    "run": solver_mode_flow,
                    "verdict": verdict_file.read_text().strip(),
                    "wallclock_seconds": float(wall_file.read_text().strip())
                    if wall_file.exists()
                    else None,
                    "exit_status": int(exit_file.read_text().strip())
                    if exit_file.exists()
                    else None,
                    "dir": str(run_dir.relative_to(baseline_dir)),
                }
            )
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("baseline_dir", nargs="?", default="tests/baseline")
    p.add_argument("--json", default=None)
    args = p.parse_args()

    baseline_dir = Path(args.baseline_dir)
    rows = collect(baseline_dir)

    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=2))

    npass = sum(1 for r in rows if r["verdict"] == "PASS")
    nfail = sum(1 for r in rows if r["verdict"] == "FAIL")
    nnocheck = sum(1 for r in rows if r["verdict"] == "NO_CHECKSUM")

    print(f"{'VERDICT':8} {'CASE':28} {'RUN':24} {'WALL(s)':>10} {'EXIT':>5}")
    for r in rows:
        wall = f"{r['wallclock_seconds']:.1f}" if r["wallclock_seconds"] is not None else "-"
        print(f"{r['verdict']:8} {r['case']:28} {r['run']:24} {wall:>10} {r['exit_status']!s:>5}")
    print()
    print(f"{npass} PASS, {nfail} FAIL, {nnocheck} NO_CHECKSUM, {len(rows)} total")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
