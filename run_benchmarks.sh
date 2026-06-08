#!/usr/bin/env bash
# run_benchmarks.sh — FloodA5 thread-scaling benchmark sweep
#
# Runs benchmark_sim.jl once per thread count, appending all results to a
# single CSV.  Equivalent to run_benchmarks.ps1 for Linux / macOS / WSL.
#
# Usage (from project root):
#   bash run_benchmarks.sh
#   bash run_benchmarks.sh --out test/my_results.csv --steps 100
#   bash run_benchmarks.sh --max-threads 8
#
# Options:
#   --out FILE          CSV output path        (default: test/benchmark_results.csv)
#   --warmup N          Warmup steps           (default: 10)
#   --steps  N          Timed steps            (default: 200)
#   --dt     S          Fixed timestep (s)     (default: 30)
#   --config FILE       benchmark_config.json  (default: built-in list)
#   --max-threads N     Override CPU detection (default: nproc)
#   --delete-existing   Delete CSV before run  (fresh results)
#   --help              Print this message

set -euo pipefail

OUT="test/benchmark_results.csv"
WARMUP=10
STEPS=200
DT=30
CONFIG=""
MAX_THREADS=0
DELETE_EXISTING=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out)             OUT="$2";          shift 2 ;;
        --warmup)          WARMUP="$2";       shift 2 ;;
        --steps)           STEPS="$2";        shift 2 ;;
        --dt)              DT="$2";           shift 2 ;;
        --config)          CONFIG="$2";       shift 2 ;;
        --max-threads)     MAX_THREADS="$2";  shift 2 ;;
        --delete-existing) DELETE_EXISTING=1; shift   ;;
        --help|-h)
            grep '^#' "$0" | sed 's/^# \{0,1\}//' | head -30
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

JULIA=$(command -v julia 2>/dev/null || true)
if [[ -z "$JULIA" ]]; then
    echo "ERROR: julia not found on PATH." >&2; exit 1
fi

if [[ "$MAX_THREADS" -le 0 ]]; then
    if command -v nproc &>/dev/null; then
        MAX_THREADS=$(nproc)
    elif [[ "$(uname)" == "Darwin" ]]; then
        MAX_THREADS=$(sysctl -n hw.logicalcpu)
    else
        MAX_THREADS=4
    fi
    echo "Detected $MAX_THREADS logical CPU cores."
fi

THREAD_COUNTS=(1)
t=2
while [[ $t -lt $MAX_THREADS ]]; do
    THREAD_COUNTS+=($t)
    t=$((t * 2))
done
if [[ "${THREAD_COUNTS[-1]}" -ne "$MAX_THREADS" ]]; then
    THREAD_COUNTS+=($MAX_THREADS)
fi

if [[ "$DELETE_EXISTING" -eq 1 && -f "$OUT" ]]; then
    rm "$OUT"; echo "Deleted existing $OUT"
fi

COMMON_ARGS=(--out "$OUT" --warmup "$WARMUP" --steps "$STEPS" --dt "$DT")
if [[ -n "$CONFIG" ]]; then COMMON_ARGS+=(--config "$CONFIG"); fi

TOTAL=${#THREAD_COUNTS[@]}
RUN_INDEX=0

echo ""; printf '%0.s=' {1..72}; echo
echo "FloodA5 benchmark sweep  —  $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Thread counts : ${THREAD_COUNTS[*]}"
echo "  Output CSV    : $OUT"
echo "  Warmup/Timed  : $WARMUP / $STEPS steps   dt=${DT}s"
printf '%0.s=' {1..72}; echo

for N in "${THREAD_COUNTS[@]}"; do
    RUN_INDEX=$((RUN_INDEX + 1))
    echo ""; echo "── Run $RUN_INDEX / $TOTAL  (--threads $N) ──"
    START=$(date +%s)
    "$JULIA" --threads "$N" benchmark_sim.jl "${COMMON_ARGS[@]}"
    EXIT_CODE=$?
    ELAPSED=$(( $(date +%s) - START ))
    if [[ $EXIT_CODE -ne 0 ]]; then
        echo "WARNING: julia exited with code $EXIT_CODE for --threads $N" >&2
    else
        echo "  Completed in ${ELAPSED}s"
    fi
done

echo ""; printf '%0.s=' {1..72}; echo
echo "Sweep complete.  Results appended to: $OUT"
printf '%0.s=' {1..72}; echo; echo ""
