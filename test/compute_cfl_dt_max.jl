#!/usr/bin/env julia
# ============================================================================
# compute_cfl_dt_max.jl
#
# Phase A, Task 4 helper (FloodA5_PhaseA_ImplementationScope.md).
#
# Prints a recommended --dt-max for a given mesh, computed from the SAME
# wave-speed CFL formula the model uses at runtime (_cfl_dt in
# surfacewater/flow2d.jl):
#
#     dt = courant * dx_min / sqrt(g * h_ref)
#
# where dx_min = sqrt(minimum(cell_area)) — the same quantity _cfl_dt itself
# computes — and h_ref is a representative reference depth (NOT measured at
# runtime; see rationale below).
#
# Why this is needed at all: dt_max is only ever an upper CEILING in the
# simulation loop (`dt = min(_cfl_dt(state, method), dt_max, ...)`), and
# _cfl_dt already recomputes a resolution-correct dt every step once real
# water depths exist. So a single fixed --dt-max would still be CFL-SAFE at
# every resolution (never a stability risk) — but it would NOT be
# CFL-CONSISTENT for the purposes of a fair resolution-convergence sweep: a
# fixed dt_max is a proportionally much looser cap at coarse resolution
# (large dx_min) than at fine resolution (small dx_min) during the dry/
# near-dry startup phase, when _cfl_dt's own fallback (60.0, when no wet
# cells exist yet) doesn't yet reflect the resolution at all. Scaling
# dt_max by dx_min, as this script does, gives every resolution in the
# sweep the same Courant-number margin at t=0, per
# FloodA5_PhaseA_ImplementationScope.md Task 4's explicit requirement.
#
# On the choice of h_ref: since this is computed BEFORE the simulation
# runs (no wet cells exist yet to measure a real depth from), h_ref is a
# single representative depth, held IDENTICAL across every resolution in
# the sweep — its absolute value doesn't need to be "correct" for
# CFL-consistency purposes (only relative dx_min scaling matters for that),
# but a wildly unrealistic h_ref would still make dt_max either uselessly
# tiny or uselessly large as a ceiling. Default 1.0 m is a middle-of-the-
# road guess for this scenario (shallow closed-basin filling test); pass
# --h-ref to override if you have a better estimate (e.g. from a prior run
# at one resolution).
#
# This calls the REAL production initialise_flow_model (StandardFlow) to
# get cell_area, rather than recomputing polygon area independently —
# guaranteeing zero risk of divergence from what the actual --meshload run
# will use (cell_area is recomputed from the polygon at init time for
# non-SGS meshes, not trusted from the parquet — see Bug 39 in
# PROJECT_STATE.md).
#
# Usage:
#   julia --project=. test/compute_cfl_dt_max.jl <mesh.parquet> [options]
#
# Options:
#   --courant N   Courant number (default 0.7, matches _cfl_dt's default)
#   --h-ref N     Reference depth in metres (default 1.0)
#   --quiet       Print only the bare numeric dt_max value (for scripting —
#                 e.g. $dtmax = julia compute_cfl_dt_max.jl mesh.parquet --quiet)
# ============================================================================

const REPO_ROOT = dirname(@__DIR__)
include(joinpath(REPO_ROOT, "mesh", "A5Grid.jl"))
include(joinpath(REPO_ROOT, "FloodModel.jl"))

function _argval(args::Vector{String}, flag::String, default=nothing)
    idx = findfirst(==(flag), args)
    idx === nothing && return default
    idx == length(args) && error("$flag requires a value")
    return args[idx + 1]
end

args_v = String.(ARGS)
QUIET  = "--quiet" in args_v

length(args_v) >= 1 && !startswith(args_v[1], "--") ||
    error("Usage: julia test/compute_cfl_dt_max.jl <mesh.parquet> [--courant N] [--h-ref N] [--quiet]")
mesh_path = args_v[1]

courant = parse(Float64, _argval(args_v, "--courant", "0.7"))
h_ref   = parse(Float64, _argval(args_v, "--h-ref", "1.0"))

mesh  = load_mesh_geoparquet(mesh_path)
state = initialise_flow_model(mesh, StandardFlow())   # cheap: no simulation loop run

dx_min = sqrt(minimum(state.cell_area))
dt_max = courant * dx_min / sqrt(9.81 * h_ref)

if QUIET
    println(round(dt_max, digits=3))
else
    println("Mesh: $mesh_path")
    println("  n_cells        = $(length(state.cell_ids))")
    println("  min cell_area  = $(round(minimum(state.cell_area), digits=1)) m²")
    println("  dx_min         = $(round(dx_min, digits=2)) m")
    println("  courant        = $courant")
    println("  h_ref          = $h_ref m  (representative depth, not measured — see script docstring)")
    println("  → recommended --dt-max = $(round(dt_max, digits=2))")
end
