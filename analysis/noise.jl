# analysis/noise.jl — checkerboard / velocity-noise metrics (WP2 diagnostic).
#
# Over "fully wet interior" cells (depth > hmin, every neighbour too, ≥ 3
# neighbours):
#   r_rms    rms of r_i = WSE_i − mean(WSE of neighbours), mm. Smooth fields
#            give O(dx²·curvature); a checkerboard gives |r| ~ its amplitude.
#   extrema  fraction of those cells that are strict local WSE extrema with
#            |r_i| > rmin (2 mm). ~0 for a smooth field, large for a
#            checkerboard (1 for a perfect two-colour one on a square grid;
#            pentagon meshes cannot be two-coloured, so less there).
#   dv_rms   rms of |v_i − mean(v of neighbours)|, m/s (velocity noise).

export noise_metrics, noise_series, noise_summary

# Neighbours of cell i: a max_nb × n matrix (0 = none) or a vector of lists.
_nbrs(nb::AbstractMatrix{<:Integer}, i) = (j for j in view(nb, :, i) if j != 0)
_nbrs(nb::AbstractVector{<:AbstractVector{<:Integer}}, i) = nb[i]
_ncells(nb::AbstractMatrix) = size(nb, 2)
_ncells(nb::AbstractVector) = length(nb)

"""
    noise_metrics(wse, h, u, v, neighbours; hmin = 0.05, rmin = 2e-3)
        -> (n, r_rms, extrema, dv_rms)

Checkerboard and velocity-noise metrics of one frame. `wse` (m), `h` depth
(m), `u`, `v` velocity components (m/s) per cell; `neighbours` either a
max_nb × n matrix (1-based, 0 = none, as `MeshGeom.neighbours`) or a vector of
neighbour index lists. Over the fully wet interior cells (`h > hmin` with a
finite WSE, ≥ 3 neighbours, all of them likewise):

- `n`: number of such cells;
- `r_rms` (mm): rms of rᵢ = WSEᵢ − mean(WSE of neighbours) — O(dx²·curvature)
  for a smooth surface, ~ the amplitude for a checkerboard;
- `extrema`: fraction that are strict local WSE extrema with |rᵢ| > `rmin` —
  ~0 smooth, large for a checkerboard (1 for a perfect ±a pattern on a
  square grid; checkerboarding over half the wet area gives ~0.5);
- `dv_rms` (m/s): rms of |vᵢ − mean(v of neighbours)|.

All NaN (and `n = 0`) if there are no such cells.
"""
function noise_metrics(wse::AbstractVector, h::AbstractVector, u::AbstractVector,
                       v::AbstractVector, nb; hmin::Real = 0.05, rmin::Real = 2e-3)
    n = length(h)
    length(wse) == length(u) == length(v) == n == _ncells(nb) ||
        error("noise_metrics: wse, h, u, v and neighbours must cover the same $n cells")
    r = Float64[]; dv = Float64[]; nextr = 0
    wn = Float64[]; un = Float64[]; vn = Float64[]
    for i in 1:n
        (h[i] > hmin && isfinite(wse[i])) || continue
        empty!(wn); empty!(un); empty!(vn)
        ok = true
        for j in _nbrs(nb, i)
            if !(h[j] > hmin && isfinite(wse[j]))
                ok = false; break
            end
            push!(wn, wse[j]); push!(un, u[j]); push!(vn, v[j])
        end
        (ok && length(wn) >= 3) || continue
        wi = wse[i]
        ri = wi - mean(wn)
        push!(r, ri)
        (abs(ri) > rmin && (all(<(wi), wn) || all(>(wi), wn))) && (nextr += 1)
        push!(dv, hypot(u[i] - mean(un), v[i] - mean(vn)))
    end
    isempty(r) && return (n = 0, r_rms = NaN, extrema = NaN, dv_rms = NaN)
    return (n = length(r), r_rms = 1e3 * sqrt(mean(r .^ 2)), extrema = nextr / length(r),
            dv_rms = sqrt(mean(dv .^ 2)))
end

"""
    noise_metrics(r::Run, k = nframes(r); sgs = (r.flow_model == "sgs"), hmin = 0.05, rmin = 2e-3)

`noise_metrics` of frame `k` of a run: WSE from `wse(r, k; sgs)`, velocities
`vel_u`, `vel_v`, the mesh's neighbours (format 1: load the run with its
mesh's `parquet`).
"""
function noise_metrics(r::Run, k::Integer = nframes(r); sgs::Bool = r.flow_model == "sgs",
                       hmin::Real = 0.05, rmin::Real = 2e-3)
    size(r.mesh.neighbours, 1) > 0 ||
        error("$(r.label): noise metrics need neighbour lists (format 1: pass the mesh's parquet)")
    sgs && !haskey(r.mesh.static, "sgs_z_min") &&
        error("$(r.label): SGS WSE needs sgs_z_min (format 1: pass the mesh's parquet)")
    h = depth(r, k)
    w = (sgs ? r.mesh.static["sgs_z_min"] : r.mesh.elevation) .+ h    # as `wse`
    return noise_metrics(w, h, field(r, "vel_u", k), field(r, "vel_v", k), r.mesh.neighbours;
                         hmin, rmin)
end

"""
    noise_series(r::Run; every = 1, sgs = (r.flow_model == "sgs"), hmin = 0.05, rmin = 2e-3)
        -> Vector{(t, k, m)}

`noise_metrics` (`m`) of every `every`-th frame `k` (k % every == 0) and the
last frame, with its time `t`.
"""
function noise_series(r::Run; every::Integer = 1, sgs::Bool = r.flow_model == "sgs",
                      hmin::Real = 0.05, rmin::Real = 2e-3)
    nf = nframes(r)
    return [(t = r.t[k], k = k, m = noise_metrics(r, k; sgs, hmin, rmin))
            for k in 1:nf if k % every == 0 || k == nf]
end

"""
    noise_summary(series) -> (last, max_r_rms, t_max_r_rms, max_extrema, max_dv_rms)

From a `noise_series`: the last entry, and over the entries with wet interior
cells (`n > 0`) the largest r_rms with its time, the largest extrema fraction
and the largest dv_rms (NaN if there are none).
"""
function noise_summary(series::AbstractVector)
    isempty(series) && error("noise_summary: empty series")
    valid = filter(s -> s.m.n > 0, series)
    imax = isempty(valid) ? 0 : argmax([s.m.r_rms for s in valid])
    return (last = series[end],
            max_r_rms = imax > 0 ? valid[imax].m.r_rms : NaN,
            t_max_r_rms = imax > 0 ? valid[imax].t : NaN,
            max_extrema = isempty(valid) ? NaN : maximum(s.m.extrema for s in valid),
            max_dv_rms = isempty(valid) ? NaN : maximum(s.m.dv_rms for s in valid))
end
