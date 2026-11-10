#!/usr/bin/env bash
# test/run_wp2_carlisle.sh
# ------------------------
# WP2 acceptance runs (FloodA5_Review_and_ForwardPlan_2026-10.md §9 WP2):
# Carlisle res-18 and res-20, standard flow, the 2005-event inflows
# (carlisle_wgs84.bci/.bdy: three QVAR inflows, free west boundary) plus
# 10 mm/hr uniform rain, 12 h, default open boundaries elsewhere.
#   A   defaults
#   B   diamond + vector friction + cell momentum            (best WP1 set)
#   C   diamond + vector + blend + positivity + local CFL    (full WP2 set)
#   C1  C with the legacy limiter          (res-18 only)
#   C2  C with the legacy CFL              (res-18 only)
# Outputs and logs go to test/wp1_runs/ (git-ignored). Analyse with
# test/wp2_checkerboard.jl. Run from the repo root of the tree to test:
#   bash test/run_wp2_carlisle.sh [parallel_jobs] [threads_per_job] [filter] [out_dir]
set -u
P="${1:-2}"; T="${2:-4}"; FILTER="${3:-.}"; OUT="${4:-test/wp1_runs}"
DATA="${FLOODA5_DATA:-test}"      # where the git-ignored meshes/BC files live
mkdir -p "$OUT"
export FLOODA5_NO_VIS=1

common="--flow-model standard --inflow-bci $DATA/carlisle/carlisle_wgs84.bci --rainfall 10 --sim-duration 43200 --dt-max 10 --output-interval 1800"
A=""
B="--gradient-correction on --face-flux-method diamond --friction-model vector --momentum-model cell"
C="--gradient-correction on --face-flux-method diamond --friction-model vector --momentum-model blend --limiter positivity --cfl-mode local"
C1="--gradient-correction on --face-flux-method diamond --friction-model vector --momentum-model blend --cfl-mode local"
C2="--gradient-correction on --face-flux-method diamond --friction-model vector --momentum-model blend --limiter positivity"

jobs=()
for res in 18 20; do
  mesh="$DATA/carlisle/carlisle_mesh${res}_standard.parquet"
  jobs+=("carlisle_res${res}_A|--meshload $mesh $common $A")
  jobs+=("carlisle_res${res}_B|--meshload $mesh $common $B")
  jobs+=("carlisle_res${res}_C|--meshload $mesh $common $C")
done
jobs+=("carlisle_res18_C1|--meshload $DATA/carlisle/carlisle_mesh18_standard.parquet $common $C1")
jobs+=("carlisle_res18_C2|--meshload $DATA/carlisle/carlisle_mesh18_standard.parquet $common $C2")

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

printf '%s\n' "${jobs[@]}" | grep -E "^[^|]*(${FILTER})" | xargs -d '\n' -P "$P" -I{} bash -c 'run_one "$@"' _ {}
