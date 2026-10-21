#!/usr/bin/env bash
# test/run_wp1_matrix.sh
# ----------------------
# WP1 test matrix (FloodA5_Review_and_ForwardPlan_2026-10.md §9 WP1):
#   Cartesian flat point spread, dx 25 / 12.5 / 6.25 m × {component, vector}
#     friction × {edge, cell} momentum; 10 m³/s central source, closed, 3600 s.
#   A5 planar slope (test/planar_embankment), res 16 / 18 × {legacy, diamond}
#     gradient × {component, vector} × {edge, cell}; 0.1 m³/s source, closed,
#     72 000 s. Plus diamond + vector + edge at θ = 1 (no Q-centring).
# Outputs and logs go to test/wp1_runs/ (git-ignored). Analyse with
# test/wp1_metrics.jl. Run from the repo root in Git Bash:
#   bash test/run_wp1_matrix.sh [parallel_jobs] [threads_per_job] [filter]
# `filter` is a grep pattern on the run name (e.g. "^cart" or "res18").
set -u
P="${1:-2}"; T="${2:-4}"; FILTER="${3:-.}"
OUT=test/wp1_runs; mkdir -p "$OUT"
export FLOODA5_NO_VIS=1

jobs=()
for dx in 25 12p5 6p25; do
  case $dx in 25) mesh=test/square/square_cartesian_dx25_flat.parquet ;;
              12p5) mesh=test/square/square_cartesian_dx12.parquet ;;
              6p25) mesh=test/square/square_cartesian_dx6p25_flat.parquet ;; esac
  for fr in component vector; do for mom in edge cell; do
    jobs+=("cart_dx${dx}_${fr}_${mom}|--meshload $mesh --flow-model standard --injection-point 0.0,0.0,10.0 --closed-boundaries --sim-duration 3600 --dt-max 2 --output-interval 600 --friction-model $fr --momentum-model $mom")
  done; done
done
for res in 16 18; do
  mesh=test/planar_embankment/planar_mesh${res}_std.parquet
  for g in legacy diamond; do
    gflags="--gradient-correction off"
    [ "$g" = diamond ] && gflags="--gradient-correction on --face-flux-method diamond"
    for fr in component vector; do for mom in edge cell; do
      jobs+=("planar_res${res}_${g}_${fr}_${mom}|--meshload $mesh --flow-model standard --injection-point 51.0001,-0.0434,0.1 --closed-boundaries --sim-duration 72000 --dt-max 2 --output-interval 300 $gflags --friction-model $fr --momentum-model $mom")
    done; done
  done
  jobs+=("planar_res${res}_diamond_vector_edge-theta1|--meshload $mesh --flow-model standard --injection-point 51.0001,-0.0434,0.1 --closed-boundaries --sim-duration 72000 --dt-max 2 --output-interval 300 --gradient-correction on --face-flux-method diamond --friction-model vector --momentum-model edge --q-centre-theta 1.0")
done

run_one() {
  name="${1%%|*}"; args="${1#*|}"
  if [ -s "$OUT/$name.h5" ] && grep -q "Simulation finished" "$OUT/$name.log" 2>/dev/null; then
    echo "skip  $name (done)"; return 0
  fi
  start=$(date +%s)
  # shellcheck disable=SC2086
  julia --project=. --threads "$T" FloodModel.jl $args --output "$OUT/$name.h5" > "$OUT/$name.log" 2>&1
  rc=$?
  echo "$( [ $rc -eq 0 ] && echo done || echo FAIL)  $name  ($(( $(date +%s) - start )) s)"
}
export -f run_one; export OUT T

printf '%s\n' "${jobs[@]}" | grep -E "^[^|]*(${FILTER})" | xargs -d '\n' -P "$P" -I{} bash -c 'run_one "$@"' _ {}
