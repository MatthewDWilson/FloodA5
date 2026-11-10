# analysis/compare.jl — run vs reference, and twin-run divergence.
#
# compare_fields / compare: one frame of a run against a reference run on the
# same mesh (from test/diag_vs_reference.jl, R1 2026-10).
# divergence_row / divergence / growth_time: how two runs that should agree
# drift apart, frame by frame (from test/diag_divergence.jl, D1 2026-10).
# The array functions keep the original scripts' arithmetic (the same
# reductions in the same order), so the wrapper scripts print identical output.

export compare_fields, compare, divergence_row, divergence, growth_time

"""
    compare_fields(h_ref, h, area; s_ref = nothing, s = nothing,
                   threshold = 0.1, wet = 0.01, deep = 3.0) -> NamedTuple

Compare depth field `h` with the reference `h_ref` (same cells, same order;
`area` in m²; `s_ref`, `s` optional cell speeds in m/s). Fields:

- `storage_ref`, `storage`: Σ h·area over cells with a finite depth (m³);
  `storage_ratio = storage / storage_ref`.
- `rms`: rms of h − h_ref (m) over the `n_wet` cells with finite depths in
  both that are wet (> `wet`, default 1 cm) in either (NaN if there are none).
- `hits`, `misses`, `false_alarms`: cells with h_ref > `threshold` and
  h > `threshold`; h_ref only; h only.
- `csi`: critical success index hits / (hits + misses + false alarms);
  `bias`: count(h > threshold) / count(h_ref > threshold). Both denominators
  are clamped to ≥ 1, so two fields dry at the threshold give 0, 0.
- `speed_ref`, `speed`: mean speed over the `n_deep` cells with h_ref ≥
  `deep` (3 m); `speed_ratio = speed / speed_ref` (NaN if no such cell or no
  speeds given).
"""
function compare_fields(h_ref::AbstractVector, h::AbstractVector, area::AbstractVector;
                        s_ref = nothing, s = nothing, threshold::Real = 0.1,
                        wet::Real = 0.01, deep::Real = 3.0)
    length(h) == length(h_ref) == length(area) ||
        error("compare_fields: lengths differ ($(length(h_ref)), $(length(h)), $(length(area)))")
    h1 = h_ref; h2 = h
    w = [i for i in eachindex(h1) if isfinite(h1[i]) && isfinite(h2[i]) && (h1[i] > wet || h2[i] > wet)]
    rms = sqrt(mean(abs2, h2[w] .- h1[w]))
    a = h1 .> threshold; b = h2 .> threshold
    hits = count(a .& b); miss = count(a .& .!b); fa = count(.!a .& b)
    d = [i for i in eachindex(h1) if h1[i] >= deep]
    st(H) = sum(H[i] * area[i] for i in eachindex(H) if isfinite(H[i]))
    sp(S) = (S === nothing || isempty(d)) ? NaN : mean(S[i] for i in d)
    st_ref = st(h1); st_run = st(h2)
    sp_ref = sp(s_ref); sp_run = sp(s)
    return (storage_ref = st_ref, storage = st_run, storage_ratio = st_run / st_ref,
            rms = rms, n_wet = length(w),
            csi = hits / max(hits + miss + fa, 1), bias = count(b) / max(count(a), 1),
            hits = hits, misses = miss, false_alarms = fa,
            speed_ref = sp_ref, speed = sp_run, speed_ratio = sp_run / sp_ref, n_deep = length(d))
end

"""
    compare(ref::Run, run::Run; k_ref = nothing, k_run = nothing,
            threshold = 0.1, area = nothing) -> NamedTuple

`compare_fields` for frame `k_run` of `run` against frame `k_ref` of `ref`,
with speeds from `speed` (|(vel_u, vel_v)|). Frames default to matched
times: if neither is given, the last time both runs reach; if one is given,
the other run's frame nearest its time. The result adds `k_ref`, `k_run`,
`t_ref` and `t` — check them: explicitly given frames need not match in time.
`area` defaults to the reference's cell areas (the run's if the reference
has none). The runs must be on the same mesh in the same cell order.
If either run is SGS, storage comes from `storage` (Σ volume); its depths
are measured from each cell's lowest point, so rms/CSI/bias against a
non-SGS run compare different quantities.
"""
function compare(ref::Run, run::Run; k_ref = nothing, k_run = nothing,
                 threshold::Real = 0.1, area = nothing)
    length(ref.mesh) == length(run.mesh) ||
        error("compare: $(ref.label) has $(length(ref.mesh)) cells, $(run.label) $(length(run.mesh))")
    if k_ref === nothing && k_run === nothing
        t = min(ref.t[end], run.t[end])
        k_ref = frame_at(ref, t); k_run = frame_at(run, t)
    elseif k_ref === nothing
        k_ref = frame_at(ref, run.t[k_run])
    elseif k_run === nothing
        k_run = frame_at(run, ref.t[k_ref])
    end
    A = area !== nothing ? area : !isempty(ref.mesh.area) ? ref.mesh.area : run.mesh.area
    isempty(A) && error("compare: no cell areas (format 1: load the runs with parquet = MESH)")
    c = compare_fields(depth(ref, k_ref), depth(run, k_run), A;
                       s_ref = speed(ref, k_ref), s = speed(run, k_run), threshold)
    if ref.flow_model == "sgs" || run.flow_model == "sgs"
        st_ref, st_run = storage(ref, k_ref), storage(run, k_run)
        c = merge(c, (storage_ref = st_ref, storage = st_run, storage_ratio = st_run / st_ref))
    end
    return merge((k_ref = k_ref, k_run = k_run, t_ref = ref.t[k_ref], t = run.t[k_run]), c)
end

"""
    divergence_row(ha, hb, ua, ub, Va, Vb; wet = 0.05, tol = 0.01) -> NamedTuple

Difference between two runs' frames (depths `h`, speeds `u`, volumes `V`)
over the cells wet (h > `wet`, default 5 cm) in either run:

- `nwet`: number of such cells;
- `rms`, `mx`: rms and maximum |ha − hb| (m), 0 if no cell is wet;
- `n1cm`: cells with |ha − hb| > `tol` (1 cm);
- `urms`: rms |ua − ub| (m/s);
- `dV`: relative difference in total stored volume (ΣVb − ΣVa) / ΣVa, sums
  over finite volumes (0 if ΣVa ≤ 0).
"""
function divergence_row(ha::AbstractVector, hb::AbstractVector, ua::AbstractVector,
                        ub::AbstractVector, va::AbstractVector, vb::AbstractVector;
                        wet::Real = 0.05, tol::Real = 0.01)
    w = findall(i -> ha[i] > wet || hb[i] > wet, eachindex(ha))
    dh = isempty(w) ? [0.0] : abs.(ha[w] .- hb[w])
    du = isempty(w) ? [0.0] : abs.(ua[w] .- ub[w])
    Va = sum(filter(isfinite, va)); Vb = sum(filter(isfinite, vb))
    return (nwet = length(w), rms = sqrt(mean(dh .^ 2)), mx = maximum(dh),
            n1cm = count(>(tol), dh), urms = sqrt(mean(du .^ 2)),
            dV = Va > 0 ? (Vb - Va) / Va : 0.0)
end

"""
    divergence(a::Run, b::Run; wet = 0.05, tol = 0.01) -> Vector{NamedTuple}

Per-frame `divergence_row` of two runs on the same mesh in the same cell
order (an error otherwise), frames paired by index (1 … the shorter run),
each row with `t` = run `a`'s frame time. Speeds: the `velocity` variable.
For twin runs (a ~1e-10 perturbation, or `--cell-order spatial` vs `file`)
and timestep-convergence pairs (the same case at different `--courant`).
"""
function divergence(a::Run, b::Run; wet::Real = 0.05, tol::Real = 0.01)
    a.mesh.ids == b.mesh.ids || error("the two files have different cell order or meshes")
    rows = NamedTuple[]
    for k in 1:min(nframes(a), nframes(b))
        r = divergence_row(field(a, "water_depth", k), field(b, "water_depth", k),
                           field(a, "velocity", k), field(b, "velocity", k),
                           field(a, "volume", k), field(b, "volume", k); wet, tol)
        push!(rows, merge((t = a.t[k],), r))
    end
    return rows
end

"""
    growth_time(rows; lo = 1e-9, hi = 1e-3) -> (tau, n)

e-folding time (s) of the rms depth difference between twin runs while it is
still small: a least-squares fit of log(rms) against t over the `divergence`
rows with `lo < rms < hi`; `tau = 1 / slope`, `n` the number of rows used
(`tau = NaN` if fewer than two).
"""
function growth_time(rows; lo::Real = 1e-9, hi::Real = 1e-3)
    sel = [r for r in rows if lo < r.rms < hi]
    length(sel) >= 2 || return (tau = NaN, n = length(sel))
    x = [r.t for r in sel]; y = log.([r.rms for r in sel])
    slope = sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^ 2)
    return (tau = 1 / slope, n = length(sel))
end
