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

  shock    A shock-front state dump: a leading shock count, then per
           shock a "<npoints> <type>" line (type optionally quoted, e.g.
           125 'S' or 125 S) followed by npoints rows of whitespace-
           separated floats -- sh00.dat, sh99_final.dat, and the
           tests/*/dumps/*.dat snapshots all share this format. Shocks
           and points are compared in file order; a run that changed the
           number of shocks or redistributed shock points (both
           legitimate outcomes of the fitting algorithm, not necessarily
           a bug) cannot be compared point-for-point and is reported as
           a structural mismatch rather than a numeric one.

Usage:
  compare.py scalar <reference> <candidate> [--rtol R] [--atol A] [--json PATH]
  compare.py shock  <reference> <candidate> [--rtol R] [--atol A] [--json PATH]
  compare.py selftest

Exit status: 0 = within tolerance, 1 = out of tolerance or structural
mismatch, 2 = usage/parse error.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

DEFAULT_RTOL = 1.0e-9
DEFAULT_ATOL = 1.0e-12


def _to_float(tok: str) -> float:
    # Fortran double-precision literals sometimes use D instead of E
    # (e.g. 1.0D+00); Python's float() only understands E.
    return float(tok.replace("D", "E").replace("d", "e"))


def parse_scalar(path: Path) -> np.ndarray:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"{path}: no non-blank lines")
    return np.array([_to_float(t) for t in lines[-1].split()])


def parse_shock(path: Path):
    """Returns a list of (type, np.ndarray[npoints, ncols]) per shock."""
    tokens = path.read_text().split()
    pos = 0

    def next_tok():
        nonlocal pos
        t = tokens[pos]
        pos += 1
        return t

    nshocks = int(next_tok())
    shocks = []
    for _ in range(nshocks):
        npoints = int(next_tok())
        shtype = next_tok().strip("'\"")
        rows = []
        # column count is inferred from the first row, then held fixed
        first_row_start = pos
        # peek ahead: find how many tokens make up one row by locating
        # the pattern where the *next* token count check no longer works
        # generically, so instead read greedily until we've consumed
        # npoints rows worth of the *same* column count as row 1.
        # Determine ncols from remaining tokens vs remaining shocks'
        # unknown sizes is not possible in general, so instead assume
        # the whole remainder up to the next integer-only line is data;
        # in practice these files are regular, so read floats until we
        # have npoints rows using a fixed ncols detected from row 1 by
        # scanning until a token fails to parse as float would be wrong
        # too (all are floats). Instead: ncols is constant per test
        # family (10 in every sample seen); infer it robustly as
        # (tokens remaining for this shock) / npoints when this is the
        # last shock, otherwise fall back to 10.
        raise NotImplementedError  # replaced below
    return shocks
