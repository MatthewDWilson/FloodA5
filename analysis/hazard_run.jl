# analysis/hazard_run.jl — hazard.jl on a loaded run (FloodA5Analysis only).

export hazard_field, statistics_maps, write_run_statistics!

"""
    hazard_field(r, var, k = nframes(r); params = HazardParams()) → Vector{Float32}

Hazard variable `var` (`HAZARD_VARS`: "depth_velocity", "hazard_rating",
"hazard_degree", "hazard_class_arr") at frame `k`, from the depth and the
cell speed (`speed`).
"""
hazard_field(r::Run, var::AbstractString, k::Integer = nframes(r); params::HazardParams = HazardParams()) =
    derive_hazard(var, depth(r, k), speed(r, k), params)

"""
    statistics_maps(r; params = HazardParams(), frames = 1:nframes(r)) → Dict{String,Vector{Float32}}

`run_statistics` over the frames `frames` of run `r` (maxima of depth, speed
and the hazard variables, time of first flooding, duration of flooding), at
the resolution of the output frames. Reads each frame once.
"""
function statistics_maps(r::Run; params::HazardParams = HazardParams(), frames = 1:nframes(r))
    ks = collect(frames)
    hv = has_field(r, "vel_u") || has_field(r, "velocity")
    return run_statistics(r.t[ks], j -> (depth(r, ks[j]), hv ? speed(r, ks[j]) : nothing),
                          length(r.mesh), params)
end

"""
    write_run_statistics!(path; flood_depth = 0.1, force = false) → Dict or nothing

Add the statistics over the run to an output file written before 0.8.13
(`/statistics`, as the model writes them at its default cadence: from the
output frames, `run_statistics` without the parameter-dependent hazard
maxima), so the viewers read them instead of deriving them on every
replay. Changes the file; files that have them already are left alone
unless `force`. Returns the statistics, or `nothing` when nothing was
written.
"""
function write_run_statistics!(path::AbstractString; flood_depth::Real = 0.1, force::Bool = false)
    FO.format_version(path) >= 2 || error("$path: format-1 files cannot take statistics (re-run them)")
    !force && FO.read_statistics(path) !== nothing && return nothing
    t = FO.frame_times(path)
    isempty(t) && error("$path has no frames")
    hv = FO.has_var(path, "velocity")
    a = h5open(path, "r") do f
        n = length(f["mesh/cell_ids"])
        acc = StatsAccumulator(n, HazardParams(flood_depth = Float64(flood_depth)); hazards = false)
        for k in eachindex(t)
            accumulate_frame!(acc, t[k], Float32.(FO.read_frame(f, "water_depth", k)),
                              hv ? Float32.(FO.read_frame(f, "velocity", k)) : nothing)
        end
        acc
    end
    st = accumulated_statistics(a)
    FO.write_statistics!(path, st, Dict{String,Any}("cadence" => "output",
        "interval_s" => length(t) > 1 ? (t[end] - t[1]) / (length(t) - 1) : 0.0,
        "flood_depth" => Float64(flood_depth), "n_states" => a.nf, "t_last" => a.t2,
        "added_by" => "FA.write_run_statistics!"))
    return st
end
