#!/bin/bash
# Runs a batch of (case solver mode flow) baseline captures in parallel,
# MAX_PARALLEL at a time. Each line of the input file (or stdin) is
# "case solver mode flow". Part of issue #11 task 4.
set -u
TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAX_PARALLEL="${MAX_PARALLEL:-5}"

running=0
while read -r case solver mode flow; do
  [ -z "$case" ] && continue
  "$TOOLS_DIR/run_baseline.sh" "$case" "$solver" "$mode" "$flow" &
  running=$((running + 1))
  if [ "$running" -ge "$MAX_PARALLEL" ]; then
    wait -n
    running=$((running - 1))
  fi
done
wait
