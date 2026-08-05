#!/bin/bash
# Runs one test case through UnDiFi-2D's fitting outer loop and archives
# the result into tests/baseline/. Part of issue #11 task 4.
#
# Usage: run_baseline.sh <case-dir-name> <solver: neo|eulfs> <mode: fitting|capturing> <flow: steady|unsteady>
set -u
CASE=$1
SOLVER=$2
MODE=$3
FLOW=$4

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CASE_DIR="$TESTS_DIR/$CASE"
OUT_DIR="$TESTS_DIR/baseline/$CASE/${SOLVER}_${MODE}_${FLOW}"
mkdir -p "$OUT_DIR"

cd "$CASE_DIR" || { echo "FAIL $CASE $SOLVER $MODE: no such test dir"; exit 1; }

# Two jobs against the SAME case directory (e.g. its fitting and capturing
# runs) must never run concurrently -- they'd race on clean.sh, log/, and
# every other file in the directory. Serialize per case-directory with an
# flock; different cases still run fully in parallel.
LOCK_FILE="$CASE_DIR/.baseline.lock"
exec 9>"$LOCK_FILE"
flock 9

# Start from a clean slate -- a leftover sh99.dat/step00501/residual_norm.dat
# etc. from a prior run in this same directory must not leak into this one.
[ -x ./clean.sh ] && ./clean.sh > /dev/null 2>&1

START=$(date +%s.%N)
bash run.sh -s "$SOLVER" -m "$MODE" -f "$FLOW" > "$OUT_DIR/stdout.log" 2>&1
STATUS=$?
END=$(date +%s.%N)
ELAPSED=$(echo "$END - $START" | bc)

echo "$ELAPSED" > "$OUT_DIR/wallclock_seconds.txt"
echo "$STATUS" > "$OUT_DIR/exit_status.txt"

# Archive whatever the run produced, keyed by solver.
CHECKSUM_FILE=""
CANDIDATE_FILE=""
if [ "$SOLVER" = "neo" ] && [ "$MODE" = "fitting" ]; then
  CHECKSUM_FILE="checksum_fitting_neo"
  [ -f residual_norm.dat ] && tail -n 1 residual_norm.dat > "$OUT_DIR/candidate_scalar.log"
  CANDIDATE_FILE="$OUT_DIR/candidate_scalar.log"
  [ -f step00501/sh99.dat ] && cp step00501/sh99.dat "$OUT_DIR/sh99.dat"
elif [ "$SOLVER" = "eulfs" ] && [ "$MODE" = "fitting" ]; then
  CHECKSUM_FILE="checksum_fitting_eulfs"
  [ -f convergenza.dat ] && tail -n 1 convergenza.dat > "$OUT_DIR/candidate_scalar.log"
  CANDIDATE_FILE="$OUT_DIR/candidate_scalar.log"
  [ -f step00501/sh99.dat ] && cp step00501/sh99.dat "$OUT_DIR/sh99.dat"
elif [ "$SOLVER" = "eulfs" ] && [ "$MODE" = "capturing" ]; then
  CHECKSUM_FILE="checksum_capturing_eulfs"
  if [ -f convhst.l2 ]; then
    tail -n 1 convhst.l2 | cut -c34- > "$OUT_DIR/candidate_scalar.log"
  fi
  CANDIDATE_FILE="$OUT_DIR/candidate_scalar.log"
fi

VERDICT="NO_CHECKSUM"
if [ "$STATUS" -ne 0 ]; then
  # A nonzero exit (now reliable: run.sh sets pipefail, so a crashed
  # binary's status survives its own `| tee run.log`) always means
  # FAIL, regardless of what a stale candidate_scalar.log left over
  # from an earlier successful run might otherwise compare as.
  VERDICT="FAIL"
elif [ -n "$CHECKSUM_FILE" ] && [ -f "$CHECKSUM_FILE" ] && [ -f "$CANDIDATE_FILE" ]; then
  if python3 "$TESTS_DIR/tools/compare.py" scalar "$CHECKSUM_FILE" "$CANDIDATE_FILE" \
      --json "$OUT_DIR/compare_report.json" > "$OUT_DIR/compare.log" 2>&1; then
    VERDICT="PASS"
  else
    VERDICT="FAIL"
  fi
fi

echo "$VERDICT" > "$OUT_DIR/verdict.txt"
printf "%-8s %-28s %-8s %-10s %-8s  wall=%6.1fs  exit=%s  %s\n" \
  "$VERDICT" "$CASE" "$SOLVER" "$MODE" "$FLOW" "$ELAPSED" "$STATUS" "$OUT_DIR"
