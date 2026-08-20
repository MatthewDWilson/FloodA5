#!/usr/bin/env julia
# ============================================================================
# count_h5_frames.jl
#
# Prints the number of frames actually present in a FloodA5 HDF5 output
# file. Needed because the number of written frames does not always equal
# the naive floor(sim_duration / output_interval) — adaptive timestepping
# and end-of-run cadence handling can produce a few fewer (observed: 228
# actual vs. 240 theoretical on a 72000s/300s run — see Task 4 sweep notes).
# Scripts that need to pick a "final" or "N% of the run" checkpoint frame
# should query this directly rather than assume the theoretical count.
#
# Usage:
#   julia --project=. test/count_h5_frames.jl <file.h5>
#
# Prints a single integer (the frame count) and nothing else.
# ============================================================================

using HDF5

length(ARGS) >= 1 || error("Usage: julia test/count_h5_frames.jl <file.h5>")
h5_path = ARGS[1]

n = h5open(h5_path, "r") do f
    length(keys(f["frames"]))
end

println(n)
