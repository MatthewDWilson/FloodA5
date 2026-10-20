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
  test/test_edge_geometry.jl
  test/test_noc_correction.jl
  test/test_sgs_unit.jl
  test/test_sgs_edge_ra_tables.jl
  test/test_inflow_point.jl
  test/test_mirror_symmetry.jl
  test/test_cell_momentum.jl
  test/test_gradient_direction.jl
  test/test_open_boundary.jl
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

echo
if [ ${#failed[@]} -eq 0 ]; then
  echo "All ${#TESTS[@]} test files passed."
else
  echo "${#failed[@]} of ${#TESTS[@]} test files FAILED:"
  printf '  %s\n' "${failed[@]}"
  exit 1
fi
