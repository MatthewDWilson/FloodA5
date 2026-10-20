# test/wp1_metrics.jl
# -------------------
# WP1 analysis of FloodA5 HDF5 output with the review's §3.4 metrics
# (FloodA5_Review_and_ForwardPlan_2026-10.md). Analysis only; reads
# /mesh/center_lons, /mesh/center_lats and /frames/*/{t,water_depth,volume}.
#
# Point spread (flat bed, central source):
#   julia --project=. test/wp1_metrics.jl spread --t 3600 --lon 0 --lat 0 FILE.h5 ...
#   A4      Σ h r² cos 4φ / Σ h r²      (< 0: diagonal elongation)
#   A6      the same with cos 6φ       (for hexagonal meshes)
#   R(h>x)  cardinal/diagonal ratio of the 8-sector maximum radius at depth > x
#
# Planar slope (north/south asymmetry about the source latitude):
#   julia --project=. test/wp1_metrics.jl planar --lat 51.0001 [--t T] FILE.h5 ...
#   asym    (V_N − V_S)/(V_N + V_S), cells within ±50 m of the source
#           latitude excluded (as test/test_planar_symmetry.jl)
#   A2      Σ h r² sin 2φ / Σ h r² about the source — signed lateral skew of
#           the plume relative to the downslope (east) axis
#   Reports the frame nearest --t (default: last), plus max|asym| over all
#   frames and the time it occurs.
#
# The frame is chosen as the one nearest the requested time, and its time is
# printed: compare runs at matched simulation times (review §3.3).

using HDF5, Printf, Statistics

const _R_EARTH = 6_371_008.8

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])

function _frames(f)
    ks = sort(collect(keys(f["frames"])))
    ts = [read(f["frames"][k]["t"]) for k in ks]
    return ks, ts
end

function _xy(f, lon0, lat0)
    lon = read(f["mesh/center_lons"]); lat = read(f["mesh/center_lats"])
    x = deg2rad.(lon .- lon0) .* _R_EARTH .* cosd(lat0)
    y = deg2rad.(lat .- lat0) .* _R_EARTH
    return x, y
end

function _harmonic(h, r, φ, k; sinus = false)
    w = h .* r .^ 2
    s = sum(w)
    s > 0 || return NaN
    return sum(w .* (sinus ? sin.(k .* φ) : cos.(k .* φ))) / s
end

function _sector_ratio(h, r, φ, thr)
    wet = findall(>(thr), h)
    length(wet) < 20 && return NaN
    sec = mod.(floor.(Int, (mod.(rad2deg.(φ[wet]), 360.0) .+ 22.5) ./ 45.0), 8)
    rad = [any(==(k), sec) ? maximum(r[wet][sec .== k]) : NaN for k in 0:7]
    return mean(rad[1:2:end]) / mean(rad[2:2:end])
end

function spread(paths; t_target, lon0, lat0)
    @printf("%-48s %8s %7s %8s %8s %9s %9s %9s\n", "file", "t (s)", "n_wet",
            "A4", "A6", "R(>1e-4)", "R(>1cm)", "R(>5cm)")
    for p in paths
        h5open(p) do f
            ks, ts = _frames(f)
            k = ks[argmin(abs.(ts .- t_target))]
            h = read(f["frames"][k]["water_depth"]); h[.!isfinite.(h)] .= 0.0
            x, y = _xy(f, lon0, lat0)
            r = hypot.(x, y); φ = atan.(y, x)
            @printf("%-48s %8.0f %7d %8.4f %8.4f %9.4f %9.4f %9.4f\n", basename(p),
                    read(f["frames"][k]["t"]), count(>(1e-4), h),
                    _harmonic(h, r, φ, 4), _harmonic(h, r, φ, 6),
                    _sector_ratio(h, r, φ, 1e-4), _sector_ratio(h, r, φ, 1e-2),
                    _sector_ratio(h, r, φ, 5e-2))
        end
    end
end

function _asym(f, k, y, x)
    v = read(f["frames"][k]["volume"]); v[.!isfinite.(v)] .= 0.0
    VN = sum(v[y .> 50.0]); VS = sum(v[y .< -50.0])
    h = read(f["frames"][k]["water_depth"]); h[.!isfinite.(h)] .= 0.0
    r = hypot.(x, y); φ = atan.(y, x)
    return (VN + VS > 1e-9 ? (VN - VS) / (VN + VS) : 0.0), _harmonic(h, r, φ, 2; sinus = true)
end

function planar(paths; t_target, lon0, lat0)
    @printf("%-52s %8s %8s %8s   %9s %8s\n", "file", "t (s)", "asym", "A2",
            "max|asym|", "at t (s)")
    for p in paths
        h5open(p) do f
            ks, ts = _frames(f)
            x, y = _xy(f, lon0, lat0)
            series = [_asym(f, k, y, x)[1] for k in ks]
            k = t_target === nothing ? ks[end] : ks[argmin(abs.(ts .- t_target))]
            a, A2 = _asym(f, k, y, x)
            imax = argmax(abs.(series))
            @printf("%-52s %8.0f %+8.4f %+8.4f   %9.4f %8.0f\n", basename(p),
                    read(f["frames"][k]["t"]), a, A2, abs(series[imax]), ts[imax])
        end
    end
end

mode  = ARGS[1]
paths = filter(a -> endswith(a, ".h5"), ARGS)
t_str = _arg("--t", nothing)
lon0  = parse(Float64, _arg("--lon", "0.0"))
lat0  = parse(Float64, _arg("--lat", "0.0"))
if mode == "spread"
    spread(paths; t_target = t_str === nothing ? Inf : parse(Float64, t_str), lon0, lat0)
elseif mode == "planar"
    planar(paths; t_target = t_str === nothing ? nothing : parse(Float64, t_str), lon0, lat0)
else
    error("mode must be 'spread' or 'planar'")
end
