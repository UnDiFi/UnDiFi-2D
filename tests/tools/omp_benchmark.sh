#!/bin/bash
# Thread-count-sweep benchmark for the OpenMP-parallelized fitting kernels
# (ROADMAP.md #16, Phase 4 exit criterion: >= 6x on 8 cores for the
# fitting-side cost). Runs the gfortran-openmp binary directly (not
# run.sh, which hardcodes the non-OpenMP gfortran-release binary and a
# fixed 501-iteration count) for a fixed, bounded iteration count -- not
# a full convergence run -- at each thread count in the list, and reports
# wall time + speedup relative to the first thread count listed (normally
# 1). A bounded count is deliberate: on refined cases (e.g.
# CircularCylinder-L3, DXCELL=0.0125) the shock-point count can approach
# the compile-time npshmax=500 cap, and the Release build has no bounds
# checking -- a bounded run gives a stable per-iteration timing average
# without that risk. See ROADMAP.md Phase 4 for the level-3 case's
# point-count-vs-iteration investigation that set this script's defaults.
#
# Usage: omp_benchmark.sh <case-dir-name> <iterations> <thread-count>[,<thread-count>...]
# Example: omp_benchmark.sh CircularCylinder-L3 100 1,2,4,8
set -u
CASE=$1
ITERS=$2
THREADS_CSV=$3

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_DIR="$(cd "$TESTS_DIR/.." && pwd)"
CASE_DIR="$TESTS_DIR/$CASE"
BIN="$REPO_DIR/build/gfortran-openmp/bin/UnDiFi-2D"
OUT_DIR="$TESTS_DIR/baseline/$CASE/omp_benchmark"
mkdir -p "$OUT_DIR"

[ -x "$BIN" ] || { echo "FAIL: $BIN not built -- run: cmake --preset gfortran-openmp && cmake --build build/gfortran-openmp"; exit 1; }
cd "$CASE_DIR" || { echo "FAIL: no such test dir $CASE_DIR"; exit 1; }

# Same per-case-directory serialization as run_baseline.sh -- this script
# and any regression run against the same case dir must not race on
# clean.sh/log/NEO_data.
LOCK_FILE="$CASE_DIR/.baseline.lock"
exec 9>"$LOCK_FILE"
flock 9

testname=$(basename "$PWD")

IFS=',' read -ra THREADS <<<"$THREADS_CSV"
declare -A WALL

for n in "${THREADS[@]}"; do
  [ -x ./clean.sh ] && ./clean.sh >/dev/null 2>&1
  cp NEO_data/textinput/inputfile-exp.txt.SF NEO_data/textinput/inputfile-exp.txt
  echo "" >NEO_data/input/vel.dat

  START=$(date +%s.%N)
  OMP_NUM_THREADS="$n" "$BIN" 0 "$ITERS" false true "$testname" >"$OUT_DIR/stdout_${n}threads.log" 2>&1
  STATUS=$?
  END=$(date +%s.%N)
  ELAPSED=$(echo "$END - $START" | bc)
  WALL[$n]=$ELAPSED

  echo "$ELAPSED" >"$OUT_DIR/wallclock_${n}threads.txt"
  printf "threads=%-3s wall=%8.2fs exit=%s\n" "$n" "$ELAPSED" "$STATUS"
done

BASE=${WALL[${THREADS[0]}]}
{
  echo
  echo "--- Summary ($CASE, $ITERS iterations, baseline=${THREADS[0]} thread(s)) ---"
  for n in "${THREADS[@]}"; do
    SPEEDUP=$(echo "scale=2; $BASE / ${WALL[$n]}" | bc)
    printf "threads=%-3s wall=%8.2fs speedup=%sx\n" "$n" "${WALL[$n]}" "$SPEEDUP"
  done
} | tee "$OUT_DIR/summary.txt"
