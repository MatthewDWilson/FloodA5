# test/diag_step_spectrum.jl
# --------------------------
# Damping diagnostics (D1, 2026-10): linear stability of the *whole* time step.
#
# The state of the standard solver is x = (V, q): cell volumes and stored face
# discharges (the cell discharge vectors used by vector friction and
# consistent Q-centring are reconstructed from q). One call of step_standard!
# is a map x → x′; its Jacobian J at a base state, by central differences of
# the real step function, contains everything — gravity operator, friction,
# limiters, Q-centring. The eigenvalues μ of J give the per-step growth
# (|μ| > 1) or damping (|μ| < 1) of each mode.
#
# Base states (closed boundaries):
#   rest   lake at rest, depth H over a flat bed (no friction to first order:
#          the pure gravity-wave problem)
#   slope  uniform depth H on a bed of slope S, after `--relax` steps with the
#          volumes held fixed, so the face discharges reach the frozen-
#          coefficient (Manning) state — friction and vector friction active
#
# Mass conservation gives one exact μ = 1 per connected closed domain; it is
# excluded from the "growing" count.
#
#   FLOODA5_NO_VIS=1 julia --project=. -t 4 test/diag_step_spectrum.jl \
#       [--mesh a5|cart] [--base rest|slope] [--depth 3.0] [--courant 0.7]
#
# D4: the last configurations add the η filter (--eta-filter C_s) to the
# positivity-limiter, θ = 1 baseline.
#
# Prints, per configuration: max |μ|, number of growing modes (|μ| > 1 + 1e-9),
# and the least-damped |μ| among the rough gravity modes: modes with at least
# 20% of their linear energy (½gAη² + ½(Lw/h)q²) in the water surface, whose
# water-surface part is in the rougher half (neighbour-difference roughness).

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using LinearAlgebra, Logging, Printf, Statistics

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const MESHKIND = _arg("--mesh", "a5")
const BASE     = _arg("--base", "rest")
const H0       = parse(Float64, _arg("--depth", "3.0"))
const COURANT  = parse(Float64, _arg("--courant", "0.7"))
const SLOPE    = parse(Float64, _arg("--slope", "1e-3"))
const RELAX    = parse(Int, _arg("--relax", "400"))

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

function _mesh()
    if MESHKIND == "cart"
        d = rad2deg(25.0 / A5Grid._EARTH_R)
        return generate_cartesian_mesh_bounds(0.0, 0.0, 22.5d, 22.5d, 25.0; lat_ref = 0.0)
    end
    aoi = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
      "coordinates":[[[172.610,-43.515],[172.630,-43.515],[172.630,-43.500],[172.610,-43.500],[172.610,-43.515]]]}}"""
    return A5Grid.mesh_for_aoi(aoi, 16)
end

function _base_state(mesh, kw)
    x = [c.center_lon for c in mesh.cells]
    z = BASE == "slope" ? SLOPE .* 111_000 .* (maximum(x) .- x) : zeros(length(x))
    mesh.static_vars["elevation"] = z                         # sills = max(bed)
    st = _quiet(() -> initialise_flow_model(mesh, StandardFlow(); kw...))
    fill!(st.ghost_cell_bc, Closed)
    st.water_depth .= H0; st.volume .= H0 .* st.cell_area
    dt = _cfl_dt(st, StandardFlow(); courant = COURANT)
    if BASE == "slope"                                         # frozen-coefficient relaxation
        V0 = copy(st.volume)
        for _ in 1:RELAX
            copyto!(st.volume, V0); fill!(st.water_depth, H0)
            step_standard!(st, dt)
        end
        copyto!(st.volume, V0); fill!(st.water_depth, H0)
    end
    _reconstruct_qvec!(st, false)
    return st, dt
end

# x → x′ through one real step
function _stepmap!(st, x, dt, n, ne)
    @views st.volume .= x[1:n]
    st.water_depth .= st.volume ./ st.cell_area                # run_simulation!'s depth sync
    @views st.edges.flux[1:ne] .= x[n+1:n+ne]
    _reconstruct_qvec!(st, false)                              # q → cell vectors (state)
    step_standard!(st, dt; velocity = false)
    return vcat(st.volume, st.edges.flux[1:ne])
end

function _jacobian(st, dt)
    n = length(st.cell_ids); ne = st.edges.n_edges
    x0 = vcat(st.volume, st.edges.flux[1:ne])
    N = n + ne; J = zeros(N, N)
    for k in 1:N
        h = 1e-6 * max(abs(x0[k]), k <= n ? st.cell_area[k] * 1e-3 : 1e-3)
        xp = copy(x0); xp[k] += h; fp = _stepmap!(st, xp, dt, n, ne)
        xm = copy(x0); xm[k] -= h; fm = _stepmap!(st, xm, dt, n, ne)
        @views J[:, k] .= (fp .- fm) ./ (2h)
    end
    return J, n, ne
end

# fraction of a mode's linear energy in the water surface
function _eta_energy(v, st, n, ne)
    η = v[1:n] ./ st.cell_area; q = v[n+1:n+ne]
    eη = sum(_G .* st.cell_area .* abs2.(η))
    eq = sum(st.edges.width[1:ne] .* st.edges.L[1:ne] ./ H0 .* abs2.(q))
    return eη / (eη + eq + eps())
end

# roughness of an eigenvector's WSE part: mean |η_i − mean(η_neighbours)| / mean|η|
function _roughness(v, st, n)
    η = abs.(v[1:n]) ./ st.cell_area
    nb = st.adj_matrix; r = 0.0; s = 0.0
    for i in 1:n
        js = [nb[k, i] for k in 1:size(nb, 1) if nb[k, i] > 0]
        isempty(js) && continue
        r += abs(η[i] - mean(η[js])); s += abs(η[i])
    end
    return s > 0 ? r / s : 0.0
end

configs = Any[
    ("recommended (consistent QC θ0.9)", scheme_options(:recommended, :standard)),
    ("recommended, θ = 1",               merge(scheme_options(:recommended, :standard), (; q_centre_theta = 1.0))),
    ("recommended, legacy QC θ0.9",      merge(scheme_options(:recommended, :standard), (; q_centre_mode = :legacy))),
    ("two-point (no GC), vector, θ = 1", merge(scheme_options(:recommended, :standard), (; gradient_correction = false, q_centre_theta = 1.0))),
    ("legacy scheme (QC θ0.9)",          scheme_options(:legacy, :standard)),
    ("legacy scheme, θ = 1",             merge(scheme_options(:legacy, :standard), (; q_centre_theta = 1.0))),
    # D4 (2026-10): the Δt-consistent baseline (positivity limiter, θ = 1) and the η filter
    ("positivity, θ = 1",                merge(scheme_options(:recommended, :standard), (; limiter = :positivity, q_centre_theta = 1.0))),
]
for cs in (0.02, 0.05, 0.08, 0.12)
    push!(configs, ("positivity, θ = 1, η filter $(cs)",
                    merge(scheme_options(:recommended, :standard), (; limiter = :positivity, q_centre_theta = 1.0, eta_filter = cs))))
end
# D6′ (2026-10): divergence-free momentum damping on top of the η filter (C_s = 0.05)
for (cm, sc) in ((0.05, :wave), (0.2, :wave), (0.5, :flow))
    push!(configs, ("pos, θ1, η 0.05, mom $(cm) $(sc)",
                    merge(scheme_options(:recommended, :standard), (; limiter = :positivity, q_centre_theta = 1.0,
                          eta_filter = 0.05, momentum_damping = cm, momentum_damping_scale = sc))))
end

mesh = _quiet(_mesh)
println("Full-step spectrum: mesh=$(MESHKIND) ($(length(mesh.cells)) cells), base=$(BASE), depth=$(H0) m, Courant=$(COURANT)",
        BASE == "slope" ? ", slope=$(SLOPE)" : "")
@printf("%-36s %8s %9s %10s %12s %14s\n", "configuration", "dt (s)", "max|μ|", "# growing", "max|μ| rough", "e-fold (steps)")
for (lab, kw) in configs
    st, dt = _base_state(deepcopy(mesh), kw)
    J, n, ne = _jacobian(st, dt)
    F = eigen(J)
    μ = abs.(F.values)
    # exclude the mass-conservation mode(s) μ = 1 (eigenvector ∝ uniform WSE)
    keep = [!(abs(F.values[k] - 1) < 1e-7 && _roughness(F.vectors[:, k], st, n) < 1e-6) for k in eachindex(μ)]
    grow = count(k -> keep[k] && μ[k] > 1 + 1e-9, eachindex(μ))
    grav = [k for k in eachindex(μ) if keep[k] && _eta_energy(F.vectors[:, k], st, n, ne) > 0.2]
    rough = Dict(k => _roughness(F.vectors[:, k], st, n) for k in grav)
    med = isempty(grav) ? 0.0 : median(collect(values(rough)))
    hi = [k for k in grav if rough[k] >= med]
    mx = maximum(μ[keep]); mxr = isempty(hi) ? NaN : maximum(μ[hi])
    ef = mx > 1 + 1e-12 ? 1 / log(mx) : Inf
    @printf("%-36s %8.3f %9.6f %10d %12.6f %14s\n", lab, dt, mx, grow, mxr,
            isfinite(ef) ? @sprintf("%.0f", ef) : "—")
end
