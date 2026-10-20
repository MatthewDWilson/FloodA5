#!/usr/bin/env bash
# test/run_jobs.sh — run FloodA5 CLI jobs from a list, N at a time.
# Each line of JOBFILE is "name|cli arguments" (blank lines and # comments
# ignored). Outputs <OUT>/<name>.h5 and <name>.log; finished jobs are skipped.
#   bash test/run_jobs.sh JOBFILE [parallel_jobs] [threads_per_job] [out_dir]
# FLOODA5_PROJECT sets the Julia project (default: .).
set -u
JOBS="$1"; P="${2:-2}"; T="${3:-4}"; OUT="${4:-test/wp1_runs}"
mkdir -p "$OUT"; export FLOODA5_NO_VIS=1
run_one() {
  name="${1%%|*}"; args="${1#*|}"
  if [ -s "$OUT/$name.h5" ] && grep -q "Simulation finished" "$OUT/$name.log" 2>/dev/null; then
    echo "skip  $name (done)"; return 0
  fi
  start=$(date +%s)
  # shellcheck disable=SC2086
  julia --project="${FLOODA5_PROJECT:-.}" --threads "$T" FloodModel.jl $args --output "$OUT/$name.h5" > "$OUT/$name.log" 2>&1
  rc=$?
  echo "$( [ $rc -eq 0 ] && echo done || echo FAIL)  $name  ($(( $(date +%s) - start )) s)"
}
export -f run_one; export OUT T
grep -v -E '^\s*(#|$)' "$JOBS" | xargs -d '\n' -P "$P" -I{} bash -c 'run_one "$@"' _ {}
