# test/channel_bench/analyse_channel.jl
# -------------------------------------
# Channel-representation series: metrics from the runs listed in a
# make_channel_meshes.jl manifest. For every case and solver, from the last
# output frame, over the measurement reach s ∈ [−0.3 L, 0.05 L] (300 m below
# the inflow, 500 m above the steep reach):
#   deta     mean water-surface error (m) in the cells that overlap the
#            channel, against the 2-D normal-flow surface Z0 − S·s + h_n.
#            Lost conveyance shows as a positive error.
#   spill    wet (> 5 cm) floodplain area per metre of reach (m) in cells
#            wholly outside the channel; 0 for in-bank flow
#   unsteady max |Δη| over the channel cells between the last two frames (m)
#   outflow  discharge leaving the domain over the last frame interval,
#            as a fraction of the inflow (1 at steady state)
# Writes RUNDIR/metrics.csv. A thin wrapper over analysis/FloodA5Analysis.jl
# (channel.jl: FA.channel_series_metrics, constants FA.CHANNEL). Format-1 runs
# take their geometry from the manifest's mesh files; format-2 runs need none.
#
#   FLOODA5_NO_VIS=1 julia --project=. test/channel_bench/analyse_channel.jl MESHDIR RUNDIR [TAG]

include(joinpath(@__DIR__, "..", "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis

const MESHDIR, RUNDIR = ARGS[1], ARGS[2]
const TAG = length(ARGS) >= 3 ? ARGS[3] : "chan"

FA.channel_series_metrics(MESHDIR, RUNDIR; tag = TAG)
