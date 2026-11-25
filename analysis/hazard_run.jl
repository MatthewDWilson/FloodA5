# analysis/hazard_run.jl — hazard.jl on a loaded run (FloodA5Analysis only).

export hazard_field, statistics_maps

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
