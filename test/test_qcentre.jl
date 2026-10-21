# test/test_qcentre.jl
# --------------------
# Q-centring on A5 (WP2c, 2026-10).
#
#   Q1  edge orientation: no Q-centring partner is reversed on a Cartesian
#       grid; on A5 (edges run from the lower to the higher cell index) a
#       large share are, and _orient_collinear! finds them
#   Q2  with oriented partners, a flow along each edge has the same sign on
#       the edge and on both partners; with legacy (unsigned) partners it
#       does not
#   Q3  frozen-coefficient conveyance on A5: θ = 1, the θ-blend and
#       consistent Q-centring reproduce the Manning friction balance without
#       deflecting the flow; legacy Q-centring (θ = 0.9, unsigned) acts as an
#       extra resistance of order (1−θ)/dt and cuts it severalfold
#   Q4  on a Cartesian grid the three Q-centring modes are the same scheme
#   Q5  consistent Q-centring keeps a lake at rest (standard and SGS kernels
#       share _q_centred_consistent; the volumetric form is checked directly)
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_qcentre.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
include(joinpath(@__DIR__, "frozen_isotropy_lib.jl"))

using Test, Logging, Statistics

const _AOI = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
  "coordinates":[[[172.610,-43.515],[172.640,-43.515],[172.640,-43.495],[172.610,-43.495],[172.610,-43.515]]]}}"""
const _A5  = A5Grid.mesh_for_aoi(_AOI, 16)
const _CART = let dlat = rad2deg(100.0 / A5Grid._EARTH_R)
    generate_cartesian_mesh_bounds(0.0, 0.0, 12dlat, 12dlat, 100.0; lat_ref = 0.0)
end
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))
_state(mesh; kw...) = _quiet(() -> initialise_flow_model(mesh, StandardFlow(); kw...))

# Fraction of partners whose stored flux has the edge's sign for a flow along
# the edge (cj → ci, i.e. against n̂_f, so q_e = +1 by the solver's convention).
function _partner_agreement(st)
    E = st.edges
    agree = 0; total = 0
    for e in 1:E.n_edges
        u = -E.nf_x[e]; v = -E.nf_y[e]                 # flow direction for edge e
        for p in (E.collinear_i[e], E.collinear_j[e])
            p == 0 && continue
            q_p = -(u * E.nf_x[abs(p)] + v * E.nf_y[abs(p)])   # same flow on the partner
            agree += flipsign(q_p, p) > 0;  total += 1
        end
    end
    return agree / total
end

@testset "Q1 — partner orientation" begin
    @test _orient_collinear!(_state(_CART).edges) == 0
    st = _state(_A5)
    n_part = count(!=(0), st.edges.collinear_i) + count(!=(0), st.edges.collinear_j)
    @test n_part > 1000
    n_rev = _orient_collinear!(st.edges)
    @test 0.3 * n_part < n_rev < 0.9 * n_part
    @test _orient_collinear!(st.edges) == 0            # idempotent
    # the init keyword does the same
    st2 = _state(_A5; q_centre_mode = :oriented)
    @test st2.edges.collinear_i == st.edges.collinear_i
    @test st2.edges.collinear_j == st.edges.collinear_j
    @test_throws ArgumentError _state(_A5; q_centre_mode = :sideways)
end

@testset "Q2 — oriented partners carry the flow's sign" begin
    a_or = _partner_agreement(_state(_A5; q_centre_mode = :oriented))
    a_lg = _partner_agreement(_state(_A5))
    @test a_or == 1.0
    @test a_lg < 0.7
    # _q_centred flips negative partners: uniform 1 on the edge and its
    # partners gives exactly 1 when they are oriented
    flux = [1.0, -1.0, 1.0]
    @test _q_centred(flux, 1, -2, 3, 0.9) ≈ 1.0
    @test _q_centred(flux, 1,  2, 3, 0.9) ≈ 0.9          # legacy reading
    @test _q_centred(flux, 1,  0, -2, 0.5) ≈ 1.0
end

@testset "Q3 — frozen-coefficient conveyance on A5" begin
    h = 0.5; n_man = 0.03
    balance = h^(5/3) / n_man                          # |q|/√S at the friction balance
    D = (; gradient_correction = true, face_flux_method = :diamond, friction_model = :vector)
    cond(; kw...) = mean(_quiet(() -> frozen_isotropy(_A5; angles_deg = [0.0, 40.0],
                                    h = h, dt = 5.0, iters = 800, kw...)).cond)
    @test cond(; D..., q_centre_theta = 1.0)          ≈ balance rtol = 0.01
    @test cond(; D..., momentum_model = :blend)       ≈ balance rtol = 0.01
    @test cond(; D..., q_centre_mode = :consistent)   ≈ balance rtol = 0.01
    @test cond(; D...) < 0.5 * balance                 # legacy Q-centring: the defect
    r = _quiet(() -> frozen_isotropy(_A5; angles_deg = [0.0, 40.0], h = h, dt = 5.0,
                                     iters = 800, D..., q_centre_mode = :consistent))
    @test r.max_abs_mean_defl < 1.0
end

# Q4/Q5 helper: a few steps on a tilted, partly wet bed
function _run(mesh; steps = 40, slope = 1e-3, depth = 0.3, kw...)
    st = _state(mesh; kw...)
    for k in eachindex(st.ghost_cell_bc); st.ghost_cell_bc[k] = Closed; end
    n = length(st.cell_ids)
    x = [c.center_lon for c in mesh.cells]
    st.elevation .= slope .* 111_000 .* (x .- minimum(x))
    for i in 1:n
        st.water_depth[i] = max(0.0, depth - st.elevation[i])
        st.volume[i] = st.water_depth[i] * st.cell_area[i]
    end
    for _ in 1:steps; step_standard!(st, 2.0); end
    return st
end

@testset "Q4 — Cartesian: legacy, oriented and consistent agree" begin
    a = _run(_CART)
    b = _run(_CART; q_centre_mode = :oriented)
    c = _run(_CART; q_centre_mode = :consistent)
    @test sum(abs, a.edges.flux) > 1e-3                 # something actually flowed
    @test b.edges.flux == a.edges.flux
    @test maximum(abs, c.edges.flux .- a.edges.flux) < 1e-9 * maximum(abs, a.edges.flux)
end

@testset "Q5 — consistent Q-centring: lake at rest, volumetric form" begin
    st = _state(_A5; q_centre_mode = :consistent)
    for k in eachindex(st.ghost_cell_bc); st.ghost_cell_bc[k] = Closed; end
    z = 2.0 .* rand(length(st.cell_ids))
    st.elevation .= z
    st.water_depth .= 5.0 .- z
    st.volume .= st.water_depth .* st.cell_area
    v0 = sum(st.volume)
    for _ in 1:50; step_standard!(st, 2.0); end
    @test maximum(abs, st.edges.flux) < 1e-10
    @test sum(st.volume) ≈ v0 rtol = 1e-12
    # volumetric (SGS) form: Q = q·w for every edge gives w_e × the unit-discharge result
    E = st.edges
    q = randn(E.n_edges); Q = q .* E.width
    qu = randn(length(st.cell_ids)); qv = randn(length(st.cell_ids))
    for e in 1:50:E.n_edges
        @test _q_centred_consistent(Q, e, E, qu, qv, 0.9; volumetric = true) ≈
              E.width[e] * _q_centred_consistent(q, e, E, qu, qv, 0.9)
    end
end

println("\nQ-centring tests complete.")
