#!/usr/bin/env python3
"""Tolerance-aware regression comparator for UnDiFi-2D outputs.

Part of issue #11 / ROADMAP.md Phase 0. Replaces byte-exact comparison
(tests/*/checksum_*, compared with `diff` in tests/checking_*_runs.sh)
with numeric L2/Linf comparison: a build, compiler, or parallelization
change that perturbs the last few digits of a converged solution should
not register as a failure, while a change that alters the physics must
reliably be caught. See tests/checking_fitting_runs.sh for the mechanism
this is meant to replace, and ROADMAP.md Part II, Phase 0 for the
"detect a deliberate 1-ulp perturbation" exit criterion this satisfies
(run `compare.py selftest` to exercise it against the fixtures already
checked into the repo).

Two file kinds are understood:

  scalar   The *last non-blank line* holds one or more whitespace-
           separated floats: checksum_fitting_neo (one value, the final
           residual_norm.dat entry), checksum_fitting_eulfs /
           checksum_capturing_eulfs (five values, the final
           convergenza.dat entry), or a live run's residual_norm.dat /
           convergenza.dat directly.

  shock    A shock-front state dump: a leading shock-count line, then
           per shock a "<npoints> <type>" line (type optionally quoted,
           e.g. 125 'S' or 125 S) followed by npoints one-line rows of
           whitespace-separated floats -- sh00.dat, sh99_final.dat, and
           the tests/*/dumps/*.dat snapshots all share this format.
           Shocks and points are compared in file order; a run that
           changed the number of shocks or the number of points on one
           (both legitimate outcomes of the fitting algorithm, not
           necessarily a bug) cannot be compared point-for-point and is
           reported as a structural mismatch rather than a numeric one.

Usage:
  compare.py scalar <reference> <candidate> [--rtol R] [--atol A] [--json PATH]
  compare.py shock  <reference> <candidate> [--rtol R] [--atol A] [--json PATH]
  compare.py selftest

Exit status: 0 = within tolerance, 1 = out of tolerance or structural
mismatch, 2 = usage/parse error.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

DEFAULT_RTOL = 1.0e-9
DEFAULT_ATOL = 1.0e-12


def _to_float(tok: str) -> float:
    # Fortran double-precision literals sometimes use D instead of E
    # (e.g. 1.0D+00); Python's float() only understands E.
    return float(tok.replace("D", "E").replace("d", "e"))


def _nonblank_lines(path: Path):
    return [ln for ln in path.read_text().splitlines() if ln.strip()]


def parse_scalar(path: Path) -> np.ndarray:
    lines = _nonblank_lines(path)
    if not lines:
        raise ValueError(f"{path}: no non-blank lines")
    return np.array([_to_float(t) for t in lines[-1].split()])


def parse_shock(path: Path):
    """Returns a list of (type, np.ndarray[npoints, ncols]) per shock,
    one entry per fitted discontinuity in the file."""
    lines = _nonblank_lines(path)
    i = 0
    nshocks = int(lines[i].split()[0])
    i += 1
    shocks = []
    for _ in range(nshocks):
        header = lines[i].split()
        i += 1
        npoints = int(header[0])
        shtype = header[1].strip("'\"") if len(header) > 1 else "?"
        rows = [
            [_to_float(t) for t in lines[i + k].split()] for k in range(npoints)
        ]
        i += npoints
        shocks.append((shtype, np.array(rows)))
    return shocks


def _within_tol(ref: np.ndarray, cand: np.ndarray, rtol: float, atol: float):
    diff = cand - ref
    l2 = float(np.linalg.norm(diff))
    linf = float(np.max(np.abs(diff))) if diff.size else 0.0
    ref_l2 = float(np.linalg.norm(ref))
    rel_l2 = l2 / ref_l2 if ref_l2 > 0 else l2
    ok = bool(np.all(np.abs(diff) <= atol + rtol * np.abs(ref)))
    return {
        "l2": l2,
        "linf": linf,
        "rel_l2": rel_l2,
        "within_tolerance": ok,
    }


def compare_scalar(ref_path: Path, cand_path: Path, rtol: float, atol: float) -> dict:
    ref = parse_scalar(ref_path)
    cand = parse_scalar(cand_path)
    if ref.shape != cand.shape:
        return {
            "kind": "scalar",
            "pass": False,
            "reason": f"column count mismatch: reference has {ref.shape[0]}, "
            f"candidate has {cand.shape[0]}",
        }
    stats = _within_tol(ref, cand, rtol, atol)
    return {
        "kind": "scalar",
        "reference": str(ref_path),
        "candidate": str(cand_path),
        "reference_values": ref.tolist(),
        "candidate_values": cand.tolist(),
        "rtol": rtol,
        "atol": atol,
        **stats,
        "pass": stats["within_tolerance"],
    }


def write_shock(path: Path, shocks) -> None:
    lines = [str(len(shocks))]
    for shtype, arr in shocks:
        lines.append(f"{arr.shape[0]} '{shtype}'")
        for row in arr:
            lines.append(" ".join(repr(float(v)) for v in row))
    path.write_text("\n".join(lines) + "\n")


def compare_shock(ref_path: Path, cand_path: Path, rtol: float, atol: float) -> dict:
    ref_shocks = parse_shock(ref_path)
    cand_shocks = parse_shock(cand_path)
    if len(ref_shocks) != len(cand_shocks):
        return {
            "kind": "shock",
            "pass": False,
            "reason": f"shock count mismatch: reference has {len(ref_shocks)}, "
            f"candidate has {len(cand_shocks)}",
        }
    per_shock = []
    overall_pass = True
    for idx, ((rtype, rarr), (ctype, carr)) in enumerate(zip(ref_shocks, cand_shocks)):
        if rarr.shape != carr.shape:
            per_shock.append(
                {
                    "shock": idx,
                    "pass": False,
                    "reason": f"point/column count mismatch: reference "
                    f"{rarr.shape}, candidate {carr.shape}",
                }
            )
            overall_pass = False
            continue
        stats = _within_tol(rarr.ravel(), carr.ravel(), rtol, atol)
        stats["shock"] = idx
        stats["type"] = rtype
        stats["npoints"] = rarr.shape[0]
        stats["ncols"] = rarr.shape[1]
        stats["pass"] = stats["within_tolerance"]
        overall_pass = overall_pass and stats["pass"]
        per_shock.append(stats)
    return {
        "kind": "shock",
        "reference": str(ref_path),
        "candidate": str(cand_path),
        "rtol": rtol,
        "atol": atol,
        "per_shock": per_shock,
        "pass": overall_pass,
    }


def _report(result: dict, json_path: str | None) -> int:
    if json_path:
        Path(json_path).write_text(json.dumps(result, indent=2))
    if result["pass"]:
        print(f"PASS  {result.get('kind')}  {result.get('candidate', '')}")
        return 0
    reason = result.get("reason")
    if reason:
        print(f"FAIL  {result.get('kind')}  {reason}")
    else:
        print(f"FAIL  {result.get('kind')}  {result.get('candidate', '')}")
        for entry in result.get("per_shock", [result]):
            if not entry.get("pass", True):
                print(
                    f"  shock {entry.get('shock', '-')}: "
                    f"l2={entry.get('l2', float('nan')):.6g} "
                    f"linf={entry.get('linf', float('nan')):.6g} "
                    f"rel_l2={entry.get('rel_l2', float('nan')):.6g} "
                    f"{entry.get('reason', '')}"
                )
    return 1


def cmd_scalar(args) -> int:
    result = compare_scalar(Path(args.reference), Path(args.candidate), args.rtol, args.atol)
    return _report(result, args.json)


def cmd_shock(args) -> int:
    result = compare_shock(Path(args.reference), Path(args.candidate), args.rtol, args.atol)
    return _report(result, args.json)


def cmd_selftest(_args) -> int:
    """Exercises both comparison paths against fixtures already checked
    into the repo: an identical copy must PASS, and a 1-ulp plus a gross
    perturbation must both FAIL. This is the concrete instance of the
    ROADMAP.md Phase 0 exit criterion."""
    import shutil
    import tempfile

    repo = Path(__file__).resolve().parents[2]
    fixtures = {
        "scalar": repo / "tests" / "CircularCylinder" / "checksum_fitting_eulfs",
        "shock": repo / "tests" / "CircularCylinder" / "su2_fitting" / "sh99_final.dat",
    }
    for kind, ref in fixtures.items():
        if not ref.exists():
            print(f"SKIP  {kind} selftest: fixture {ref} not found")
            continue
        compare = compare_scalar if kind == "scalar" else compare_shock

        with tempfile.TemporaryDirectory() as td:
            td = Path(td)

            identical = td / "identical"
            shutil.copy(ref, identical)
            r = compare(ref, identical, DEFAULT_RTOL, DEFAULT_ATOL)
            assert r["pass"], f"{kind}: identical copy must PASS, got {r}"
            print(f"OK    {kind}: identical copy correctly PASSes")

            # Perturb through the parser/serializer round-trip so this
            # targets the last *parsed* value regardless of any trailing
            # content the format ignores (e.g. sh99_final.dat carries a
            # special-point section after the shock data this comparator
            # doesn't interpret).
            if kind == "scalar":
                ref_vals = parse_scalar(ref)

                def make(delta):
                    v = ref_vals.copy()
                    v[-1] = v[-1] + delta if delta != "ulp" else math.nextafter(v[-1], v[-1] + 1)
                    return v

                def write(path, vals):
                    path.write_text(" ".join(repr(float(x)) for x in vals) + "\n")

            else:
                ref_shocks = parse_shock(ref)

                def make(delta):
                    shocks = [(t, a.copy()) for t, a in ref_shocks]
                    last_type, last_arr = shocks[-1]
                    v = last_arr[-1, -1]
                    last_arr[-1, -1] = math.nextafter(v, v + 1) if delta == "ulp" else v + delta
                    return shocks

                def write(path, shocks):
                    write_shock(path, shocks)

            ulp = td / "ulp"
            write(ulp, make("ulp"))
            r = compare(ref, ulp, rtol=0.0, atol=0.0)
            assert not r["pass"], f"{kind}: 1-ulp perturbation with rtol=atol=0 must FAIL, got {r}"
            print(f"OK    {kind}: 1-ulp perturbation correctly FAILs at rtol=atol=0")
            r = compare(ref, ulp, DEFAULT_RTOL, DEFAULT_ATOL)
            assert r["pass"], f"{kind}: 1-ulp perturbation must PASS at default tolerance, got {r}"
            print(f"OK    {kind}: 1-ulp perturbation correctly PASSes at default tolerance")

            gross = td / "gross"
            write(gross, make(1.0))
            r = compare(ref, gross, DEFAULT_RTOL, DEFAULT_ATOL)
            assert not r["pass"], f"{kind}: gross perturbation must FAIL, got {r}"
            print(f"OK    {kind}: gross (+1.0) perturbation correctly FAILs")

    print("selftest: all checks passed")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    for name, fn in (("scalar", cmd_scalar), ("shock", cmd_shock)):
        sp = sub.add_parser(name)
        sp.add_argument("reference", type=str)
        sp.add_argument("candidate", type=str)
        sp.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
        sp.add_argument("--atol", type=float, default=DEFAULT_ATOL)
        sp.add_argument("--json", type=str, default=None)
        sp.set_defaults(func=fn)

    sp = sub.add_parser("selftest")
    sp.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    try:
        return args.func(args)
    except (ValueError, IndexError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
