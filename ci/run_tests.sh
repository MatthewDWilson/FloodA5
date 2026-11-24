#!/usr/bin/env bash
# Run the FloodA5 CI test set headlessly. Usage (from repo root):
#   bash ci/run_tests.sh [julia-project-dir]      # default: ci
# Runs every test file even if one fails, then reports and exits non-zero
# if any failed.
set -u
PROJECT="${1:-ci}"
export FLOODA5_NO_VIS=1
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-2}"

TESTS=(
  test/test_wp0_regressions.jl
  test/test_pya5_compat.jl
  test/test_mesh_families.jl
  test/test_dggs_meshes.jl
  test/test_dem_registration.jl
  test/test_sgs_general.jl
  test/test_default_scheme.jl
  test/test_qcentre.jl
  test/test_wp4_equivalence.jl
  test/test_wp4_stage2.jl
  test/test_limiter_stats.jl
  test/test_eta_filter.jl
  test/test_momentum_damping.jl
  test/test_swe_reference.jl
  test/test_backend.jl
  test/test_precision.jl
  test/test_wp1_vector_friction.jl
  test/test_wp2_momentum_limiters.jl
  test/test_edge_geometry.jl
  test/test_noc_correction.jl
  test/test_sgs_unit.jl
  test/test_sgs_edge_ra_tables.jl
  test/test_sgs_connectivity.jl
  test/test_inflow_point.jl
  test/test_mirror_symmetry.jl
  test/test_cell_momentum.jl
  test/test_gradient_direction.jl
  test/test_open_boundary.jl
  test/test_output_format.jl
  test/test_analysis.jl
  test/test_replay.jl
  test/test_vis_server.jl
  test/test_vis_live.jl
  test/test_plots.jl
  test/test_cases.jl
  test/test_notebooks.jl
)

failed=()
for t in "${TESTS[@]}"; do
  echo "::group::$t"
  start=$(date +%s)
  if julia --project="$PROJECT" --color=yes "$t"; then
    echo "PASS  $t  ($(( $(date +%s) - start ))s)"
  else
    echo "::error file=$t::FAILED"
    failed+=("$t")
  fi
  echo "::endgroup::"
done

# The web viewer's JavaScript (protocol, geometry, colour maps): Node's own
# test runner, no npm packages. Skipped where Node is not installed.
total=${#TESTS[@]}
if command -v node >/dev/null 2>&1; then
  total=$((total + 1))
  echo "::group::visualisation/cesium (node --test)"
  if (cd visualisation/cesium && node --test test/*.test.js); then   # file list: Node 20 takes no directory
    echo "PASS  visualisation/cesium/test"
  else
    echo "::error file=visualisation/cesium/test::FAILED"
    failed+=("visualisation/cesium/test (node)")
  fi
  echo "::endgroup::"
else
  echo "SKIP  visualisation/cesium/test (node not found)"
fi

echo
if [ ${#failed[@]} -eq 0 ]; then
  echo "All $total test files passed."
else
  echo "${#failed[@]} of $total test files FAILED:"
  printf '  %s\n' "${failed[@]}"
  exit 1
fi
