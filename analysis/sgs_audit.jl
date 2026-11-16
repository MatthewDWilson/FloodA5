# analysis/sgs_audit.jl — SGS table audit: hypsometry, sub-basins, Courant.
#
# Pure functions of arrays (no mesh loader, no DEM reader); the CLI that feeds
# them is test/diag_sgs_tables_audit.jl. Plan: projectinstructions/
# FloodA5_SGS_Conveyance_Fix_Plan_2026-11.md, Phase 0 items 2–4.
#
#   HypsCurve                    exact level-pool curve of a set of equal-area terrain samples
#   hyps_volume / hyps_wet_area / hyps_wse
#   uniform_level_table          the table build_sgs_tables! makes from sorted samples
#   sgs_table_volume / sgs_table_wse / sgs_table_wet_area
#                                the stored table as the model reads it
#   basin_merge_tree             sub-basins of a gridded cell (sublevel-set merge tree)
#   significant_basins, lowest_point_basin
#   courant_exposure             Courant numbers at a given step, and the step a depth requires
#   dist_summary                 n, median, p95, max … of a vector

export HypsCurve, hyps_volume, hyps_wet_area, hyps_wse, uniform_level_table,
       sgs_table_volume, sgs_table_wse, sgs_table_wet_area,
       Basin, BasinTree, basin_merge_tree, significant_basins, lowest_point_basin,
       courant_exposure, dist_summary

# ───────────────────────── reference hypsometry ────────────────────────────

"""
    HypsCurve(z, cell_area)

Level-pool storage curve of terrain samples `z` (any order; non-finite values
dropped) that each stand for `cell_area / n` of plan area. Exact for those
samples: V(η) = dA·Σ max(η − z, 0), with no level discretisation.
"""
struct HypsCurve
    z::Vector{Float64}       # sorted
    csum::Vector{Float64}    # prefix sums of z
    vz::Vector{Float64}      # volume at η = z[k]
    dA::Float64
    area::Float64
end

function HypsCurve(z::AbstractVector{<:Real}, cell_area::Real)
    zs = sort!(Float64[v for v in z if isfinite(v)])
    isempty(zs) && error("HypsCurve: no finite samples")
    dA = cell_area / length(zs)
    cs = cumsum(zs)
    vz = [dA * (k * zs[k] - cs[k]) for k in eachindex(zs)]
    return HypsCurve(zs, cs, vz, dA, Float64(cell_area))
end

"Stored volume (m³) below level η."
function hyps_volume(c::HypsCurve, η::Real)
    k = searchsortedlast(c.z, η)
    k == 0 && return 0.0
    return c.dA * (k * η - c.csum[k])
end

"Wetted plan area (m²) at level η (samples at or below η, as the table build counts them)."
hyps_wet_area(c::HypsCurve, η::Real) = c.dA * searchsortedlast(c.z, η)

"Level (m) at which the curve holds volume V (the lowest sample for V ≤ 0)."
function hyps_wse(c::HypsCurve, V::Real)
    V <= 0 && return c.z[1]
    k = searchsortedlast(c.vz, V)          # vz is non-decreasing; vz[1] = 0
    return c.z[k] + (V - c.vz[k]) / (c.dA * k)
end

# ───────────────────────── the model's table ───────────────────────────────

"""
    uniform_level_table(z_sorted, cell_area; n_bins = 100) -> (elev_bins, vol_curve, area_curve)

The cell table `A5Grid.build_sgs_tables!` builds from sorted, finite sample
elevations: `n_bins` levels uniformly spaced from the lowest to the highest
sample, cumulative volume and wetted area at each (same operations in the
same order, so a rebuild from the same samples reproduces the stored table).
"""
function uniform_level_table(z_sorted::AbstractVector{<:Real}, cell_area::Real; n_bins::Integer = 100)
    n = length(z_sorted)
    n > 0 || error("uniform_level_table: no samples")
    bins = range(Float64(z_sorted[1]), Float64(z_sorted[end]), length = n_bins)
    dA = cell_area / n
    eb = Vector{Float64}(undef, n_bins); vc = similar(eb); ac = similar(eb)
    cum_area = 0.0; cum_vol = 0.0; idx = 1; n_wet = 0; prev = Float64(z_sorted[1])
    for k in 1:n_bins
        w = Float64(bins[k])
        cum_vol += n_wet * dA * (w - prev)
        while idx <= n && z_sorted[idx] <= w
            cum_area += dA
            cum_vol += dA * (w - z_sorted[idx])
            n_wet += 1
            idx += 1
        end
        eb[k] = w; ac[k] = cum_area; vc[k] = cum_vol
        prev = w
    end
    return eb, vc, ac
end

"""
    sgs_table_volume(elev_bins, vol_curve, cell_area, η)

Volume the stored table holds at level η, i.e. the inverse of the model's
`wse_from_volume`: linear between levels, the full cell area above the top one.
"""
function sgs_table_volume(eb::AbstractVector, vc::AbstractVector, cell_area::Real, η::Real)
    η <= eb[1] && return 0.0
    η >= eb[end] && return vc[end] + (η - eb[end]) * cell_area
    k = searchsortedlast(eb, η)
    return vc[k] + (η - eb[k]) / (eb[k+1] - eb[k] + 1e-15) * (vc[k+1] - vc[k])
end

"Level for volume V, as `A5Grid.wse_from_volume` (same arithmetic)."
function sgs_table_wse(eb::AbstractVector, vc::AbstractVector, cell_area::Real, V::Real)
    V <= 0.0 && return Float64(eb[1])
    V >= vc[end] && return eb[end] + (V - vc[end]) / cell_area
    lo, hi = 1, length(vc)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        vc[mid] <= V ? (lo = mid) : (hi = mid)
    end
    frac = (V - vc[lo]) / (vc[hi] - vc[lo] + 1e-15)
    return eb[lo] + frac * (eb[hi] - eb[lo])
end

"Wetted plan area at level η, as `A5Grid.wetted_area_from_wse`."
function sgs_table_wet_area(eb::AbstractVector, ac::AbstractVector, cell_area::Real, η::Real)
    η <= eb[1] && return 0.0
    η >= eb[end] && return Float64(cell_area)
    k = searchsortedlast(eb, η)
    return ac[k] + (η - eb[k]) / (eb[k+1] - eb[k] + 1e-15) * (ac[k+1] - ac[k])
end

# ───────────────────────────── sub-basins ──────────────────────────────────

"""
    Basin

One branch of a cell's merge tree: a connected region of ground below a
rising level, from its own lowest point until it meets a larger region.

- `z_min`: lowest point of the region when it merged (m)
- `merge_level`: level at which it met a larger region (m); `Inf` for the root
- `area_frac`: its plan area at that level / the cell's (1 for the root)
- `depth`: `merge_level − z_min` (the root: highest point − lowest)
- `parent`: index of the region it merged into (0 for the root)
- `main_depth`, `main_area_frac`: the same two measures of the larger region
  at that merge (NaN for the root)
"""
struct Basin
    z_min::Float64
    merge_level::Float64
    area_frac::Float64
    depth::Float64
    parent::Int
    main_depth::Float64
    main_area_frac::Float64
end

"`basins`, the index `root` of the main basin, and `n`, the number of grid points used."
struct BasinTree
    basins::Vector{Basin}
    root::Int
    n::Int
end

"""
    basin_merge_tree(z, mask = isfinite.(z); flat_tol = 1e-3) -> BasinTree

Sub-basins of the gridded terrain `z` (a matrix on a regular grid; only points
with `mask` true and a finite `z` are used). As the level η rises, the ground
below η forms 4-connected regions (so a wall one grid point thick along a
diagonal holds). A region starts at a local minimum; when two regions meet,
the one with the smaller plan area at that level ends there and is recorded
with its area and depth, and the larger continues. The region that survives
every merge is the **main basin** (`root`): it is the larger side at each of
its merges and ends as the whole cell. Basins whose `parent` is the root
merged directly with the main basin; anything that had merged into them
earlier is counted in their area and depth.

Equal elevations are taken in grid order, so a flat produces zero-depth
branches; filter by depth (see `significant_basins`). So that a plateau filled
in grid order cannot pass for the main basin, a region shallower than
`flat_tol` (m) at a merge ends there whatever its area if the other has depth.
"""
function basin_merge_tree(z::AbstractMatrix{<:Real}, mask::AbstractMatrix{Bool} = isfinite.(z);
                          flat_tol::Real = 1e-3)
    nx, ny = size(z)
    pts = [i for i in LinearIndices(z) if mask[i] && isfinite(z[i])]
    n = length(pts)
    n > 0 || error("basin_merge_tree: no grid points")
    sort!(pts; by = i -> z[i], alg = MergeSort)       # stable: ties in grid order
    uf = zeros(Int, nx * ny)                          # 0 = not yet below the level
    count = zeros(Int, nx * ny)                       # at set representatives
    zlow = fill(Inf, nx * ny)
    branch = zeros(Int, nx * ny)
    find(i) = begin
        while uf[i] != i
            uf[i] = uf[uf[i]]; i = uf[i]
        end
        i
    end
    bz = Float64[]; bm = Float64[]; ba = Float64[]; bd = Float64[]; bp = Int[]
    bmd = Float64[]; bma = Float64[]
    roots = Int[]
    for p in pts
        zp = Float64(z[p])
        ix = (p - 1) % nx + 1; iy = (p - 1) ÷ nx + 1
        empty!(roots)
        for (dx, dy) in ((-1, 0), (1, 0), (0, -1), (0, 1))
            jx = ix + dx; jy = iy + dy
            (1 <= jx <= nx && 1 <= jy <= ny) || continue
            q = (jy - 1) * nx + jx
            uf[q] == 0 && continue
            r = find(q)
            r in roots || push!(roots, r)
        end
        if isempty(roots)
            uf[p] = p; count[p] = 1; zlow[p] = zp
            push!(bz, zp); push!(bm, Inf); push!(ba, 1.0); push!(bd, 0.0); push!(bp, 0)
            push!(bmd, NaN); push!(bma, NaN)
            branch[p] = length(bz)
            continue
        end
        # survivor: the largest region; a region still flat (a plateau being
        # filled in grid order, not a pool) never outlasts one with depth
        s = roots[1]
        for r in roots
            ds = zp - zlow[s] > flat_tol; dr = zp - zlow[r] > flat_tol
            ((dr && !ds) || (dr == ds && count[r] > count[s])) && (s = r)
        end
        for r in roots
            r == s && continue
            b = branch[r]
            bz[b] = zlow[r]; bm[b] = zp; ba[b] = count[r] / n; bd[b] = zp - zlow[r]
            bp[b] = branch[s]; bmd[b] = zp - zlow[s]; bma[b] = count[s] / n
            uf[r] = s; count[s] += count[r]; zlow[s] = min(zlow[s], zlow[r])
        end
        uf[p] = s; count[s] += 1
    end
    r = find(pts[1])                                   # unmerged regions (several only if the
    for p in pts                                       # mask is disconnected): the largest is the root
        uf[p] == p || continue
        bz[branch[p]] = zlow[p]
        count[p] > count[r] && (r = p)
    end
    b = branch[r]
    bd[b] = Float64(z[pts[end]]) - zlow[r]
    return BasinTree([Basin(bz[i], bm[i], ba[i], bd[i], bp[i], bmd[i], bma[i]) for i in eachindex(bz)], b, n)
end

"""
    significant_basins(tree; area, depth, two_sided = false) -> Int

Number of basins that merged directly with the main basin holding more than
`area` of the cell's plan area at that level and deeper than `depth` (m) below
it. `two_sided = true` also requires the main basin to be deeper than `depth`
below the merge level (two pools, not a pool beside a shelf).
"""
function significant_basins(t::BasinTree; area::Real, depth::Real, two_sided::Bool = false)
    return count(b -> b.parent == t.root && b.area_frac > area && b.depth > depth &&
                      (!two_sided || b.main_depth > depth), t.basins)
end

"""
    lowest_point_basin(tree) -> nothing or Basin

`nothing` if the cell's lowest point lies in the main basin; otherwise the
basin, among those that merged directly with the main one, that holds it.
"""
function lowest_point_basin(t::BasinTree)
    b = 1                                              # the first branch starts at the lowest point
    b == t.root && return nothing
    while t.basins[b].parent != t.root
        t.basins[b].parent == 0 && return t.basins[b]  # a disconnected piece of the mask
        b = t.basins[b].parent
    end
    return t.basins[b]
end

# ───────────────────────────── time step ───────────────────────────────────

"""
    courant_exposure(h, ℓ, dt; courant = 0.7, g = 9.81, wet = 1e-4)

Local Courant numbers C_i = dt·√(g h_i)/ℓ_i of the cells with h_i > `wet`
(the form of the model's `_cfl_dt_local`), for depth measure `h` (m), cell
lengths `ℓ` (m) and a step `dt` (s). Returns `(C, n, median, p95, p99, max,
share_above_0p7, share_above_1, dt_required)`; `dt_required` = courant / max
√(g h)/ℓ, the step the local rule would take with this depth.
"""
function courant_exposure(h::AbstractVector, ℓ::AbstractVector, dt::Real;
                          courant::Real = 0.7, g::Real = 9.81, wet::Real = 1e-4)
    C = [dt * sqrt(g * h[i]) / ℓ[i] for i in eachindex(h) if isfinite(h[i]) && h[i] > wet && ℓ[i] > 0]
    isempty(C) && return (C = C, n = 0, median = NaN, p95 = NaN, p99 = NaN, max = NaN,
                          share_above_0p7 = NaN, share_above_1 = NaN, dt_required = NaN)
    return (C = C, n = length(C), median = median(C), p95 = quantile(C, 0.95), p99 = quantile(C, 0.99),
            max = maximum(C), share_above_0p7 = count(>(0.7), C) / length(C),
            share_above_1 = count(>(1.0), C) / length(C), dt_required = courant * dt / maximum(C))
end

"n, min, median, mean, p95, p99 and max of the finite values of `v` (a Dict; n = 0 if none)."
function dist_summary(v::AbstractVector)
    w = Float64[x for x in v if isfinite(x)]
    isempty(w) && return Dict{String,Any}("n" => 0)
    return Dict{String,Any}("n" => length(w), "min" => minimum(w), "median" => median(w), "mean" => mean(w),
                            "p95" => quantile(w, 0.95), "p99" => quantile(w, 0.99), "max" => maximum(w))
end
