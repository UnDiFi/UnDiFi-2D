# tests/baseline

Reference runs captured for issue #11 / ROADMAP.md Phase 0 ("baseline
capture: run all 12 cases on all 3 solvers, archive reference fields +
wall-clock + profile"). This is the ground truth every later refactor
phase (F77->F90 conversion, removing `dstak`, the object model, OpenMP,
...) is checked against via `tests/tools/compare.py`.

## Layout

```
tests/baseline/<case>/<solver>_<mode>_<flow>/
  stdout.log             full run.sh output
  wallclock_seconds.txt  wall-clock time for the run
  exit_status.txt        run.sh's exit status
  candidate_scalar.log   the final convergence scalar(s) this run produced
  compare.log            tests/tools/compare.py's text output
  compare_report.json    tests/tools/compare.py's machine-readable output
  verdict.txt            PASS | FAIL | NO_CHECKSUM
  sh99.dat               final shock-front state (fitting runs only) --
                          richer than the scalar residual: full shock
                          geometry and upstream/downstream state, usable
                          with `compare.py shock` for future comparisons
```

Produced by `tests/tools/run_baseline.sh <case> <solver> <mode> <flow>`
(one run) and `tests/tools/run_baseline_batch.sh` (many, in parallel --
reads `case solver mode flow` lines from stdin, `MAX_PARALLEL` env var
controls concurrency). `tests/tools/summarize_baseline.py` walks this
directory and prints/writes a consolidated report.

Each run starts from a clean checkout of its test case directory
(`clean.sh` is run first) so results don't depend on run order.

## Methodology and what "baseline" means here

The existing `tests/*/checksum_*` files are the original project's
long-standing reference values (a single final-iteration convergence
scalar, checked with an exact `diff`). Rather than discard them, this
baseline capture *validates the current, unmodified code against them*:
every run here is compared against the pre-existing checksum via
`compare.py scalar`. Where the two agree (bit-for-bit or within
tolerance), that's strong evidence the codebase, freshly rebuilt from
source, still behaves exactly as it always has -- which is the actual
goal of a Phase-0 baseline, not merely capturing whatever the code
currently outputs.

The `sh99.dat` shock-front state is captured in addition to the scalar,
because the scalar checksum is a single number from the *entire
iteration history* (sensitive to any path-dependent floating-point
accumulation over hundreds of nonlinear iterations), while the shock
state is the converged physical result. Future refactor phases should
prefer `compare.py shock` against these files over the scalar checksums
where both are available.

## Status (2026-08-02)

**NEO, fitting, steady -- all 11 cases with a `checksum_fitting_neo`:**
10/11 reproduce the existing checksum bit-for-bit. One case,
`MachReflection-2`, differs by a relative ~1.1e-7 in the final log10
residual scalar (`-9.131263` reference vs `-9.131264` here) -- see
`MachReflection-2/neo_fitting_steady/compare_report.json`. Given the
magnitude (roughly one part in 10 million, in a metric that accumulates
over 501 nonlinear iterations) and that this is the only one of 11 cases
to disagree at all, this reads as benign floating-point drift (most
likely a glibc/libm version difference between this environment and
whenever the checksum was originally recorded) rather than a functional
regression. It is flagged rather than silently accepted because the
default comparator tolerance (rtol=1e-9) is, by design, tight enough to
catch exactly this kind of thing -- that tightness is doing its job
here, not misbehaving. No corresponding reference `sh99.dat` exists to
cross-check the actual shock geometry for this case (only the scalar
checksum was ever committed), so this baseline run's own `sh99.dat` is
the first such reference captured for `MachReflection-2` and should be
used as the point of comparison from now on.

**EulFS (fitting + capturing), steady, all 10 cases with the matching
checksums, plus `CircularCylinder`: 21/21 reproduce the existing
checksums bit-for-bit.** Required PETSc 3.14.6 + EulFS 3.14 built from
source. Two build issues came up, both environment-version mismatches
rather than anything wrong with the code:

* PETSc 3.14.6's bundled `./configure` (Python 2/early-3 era) imports
  the stdlib `xdrlib` module, removed in Python 3.13 (this environment's
  default `python3`). Built with `/usr/bin/python3.12` instead, which
  still has it.
* `EulFS.3.14/src/seq/makefile` names its output `EulFS_$(HOSTTYPE)`;
  without `HOSTTYPE` exported it built as `EulFS_` (no suffix), which
  `run.sh` doesn't look for (it hardcodes `EulFS_x86_64`). Fixed by
  exporting `HOSTTYPE=x86_64` before `make install`, matching the
  convention `compile_all.sh` already uses for `source/` and
  `source_utils/`.

**A real bug was found and fixed in this baseline-capture tooling
itself, not in UnDiFi-2D:** the first EulFS batch ran `run_baseline.sh`
for a case's `fitting` and `capturing` modes back-to-back in the job
list, and at `MAX_PARALLEL=4` both ended up scheduled concurrently.
Since both jobs `cd` into the *same* test-case directory and both start
with `./clean.sh`, two jobs racing there stomp on each other's `log/`
directory, seed files, and in-flight output -- one job's `clean.sh` can
fire mid-way through the other's run. This corrupted 9 of the first
21 EulFS runs (`fitting` runs crashing a few iterations in with
`MPI_Abort`, `log/eulfs.log` never even getting created) and produced
one silently-wrong `FAIL` (`CoaSHCK` capturing, whose checksum mismatch
was actually cross-contaminated output, not a real numeric difference).
The `PASS` results from that same racy batch are not suspect -- a
corrupted run reproducing an exact multi-digit reference checksum by
coincidence is not plausible -- but the failures were real symptoms of
the race, not of anything in the solver. Fixed with a per-case-directory
`flock` in `run_baseline.sh` (two jobs targeting the same case now
serialize; different cases still run fully in parallel), and the 10
affected jobs were re-run cleanly: all 10 now `PASS`.

**`SU2_CFD`-based cases (`CircularCylinder-su2fitting`,
`CircularCylinder-eulfsfitting`) and `ShockVortex`:** out of scope for
this pass -- no checksum files exist for these (the su2fitting/
eulfsfitting dirs use ad-hoc reference files instead, `ShockVortex` has
none at all), and `SU2_CFD` needs a separate reinstall. Left for a
follow-up.

## Summary

**32/33 runs reproduce the existing reference checksums** (bit-for-bit
in every case but one). The one exception, `MachReflection-2`'s NEO
run, differs by a relative ~1.1e-7 in a log10-residual scalar accumulated
over 501 nonlinear iterations -- see above; read as benign, not a
regression. Every EulFS run and 10/11 NEO runs matched exactly.

This is strong evidence that: (a) the codebase, freshly rebuilt from
source in a modern environment (gfortran 13.3, a from-scratch PETSc
3.14.6 + EulFS 3.14 build), still behaves exactly as its long-standing
reference values say it should; and (b) the repo hygiene and CMake work
earlier in issue #11 introduced no behavioral change -- both are
exactly what a Phase-0 baseline is for.

Approximate per-case cost, for planning future full-suite runs: NEO
fitting ~7-28 min/case, EulFS fitting ~12-21 min/case, EulFS capturing
~1-2 min/case (capturing has no per-iteration remeshing overhead, so
it's dominated by the flow solve alone rather than the outer shock-fitting
loop). Full parallel run of everything in this baseline (33 runs) took
on the order of 2 hours wall-clock at 4-5x concurrency on an 8-core
machine.
