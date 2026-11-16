# test/test_wp4_equivalence.jl
# ----------------------------
# WP4 stage 1 (2026-10): every fast path added for performance must give
# exactly (===) the result of the code it replaced.
#
#   P1  edge_area_radius == (flow_area_from_wse, hydraulic_radius_from_wse),
#       including the dry, overfull and no-table cases
#   P2  per-vertex diamond values (diamond_vertex_eta!/_wet!, diamond_dWSE_n_v,
#       diamond_stencil_wet_v) == the per-edge Dict lookups, on an A5 mesh
#   P3  step with velocity = false, then _compute_velocity!, == velocity = true;
#       the step workspace is reused (no per-step growth)
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_wp4_equivalence.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging, Random

Random.seed!(20261022)

@testset "P1 — face_apk: one lookup, exact at its levels, continuous, dry slots" begin
    # (Until 0.7.0 P1 pinned edge_area_radius to the separate per-cell lookups
    # it replaced; those are gone with the per-cell edge curves.)
    nf = A5Grid.SGS_FACE_LEVELS
    z0 = fill(NaN, 6); zr = zeros(6); A = zeros(nf, 6); P = zeros(nf, 6); K = zeros(nf, 6)
    profs = [sort(rand(40)) .* 3.0 .+ s for s in 1:5]      # slot 6 has no table
    for s in 1:5
        z0[s], zr[s], A[:, s], P[:, s], K[:, s] = A5Grid.face_table_from_profile(profs[s], 0.5)
    end
    t = SGSTable([0.0, 1.0], [0.0, 1.0], [0.0, 1.0], 100.0, 0.0, 1.0, z0, zr, A, P, K)
    for s in 1:5
        for k in 1:nf                                       # the stored value at each level
            η = k == nf ? z0[s] + zr[s] : z0[s] + zr[s] * ((k - 1) / (nf - 1))^2
            a, p, c = A5Grid.face_apk(t, s, η)
            @test a ≈ A[k, s] atol = 1e-9 * (1 + A[nf, s])
            @test c ≈ K[k, s] atol = 1e-9 * (1 + K[nf, s])
        end
        top = z0[s] + zr[s]                                 # continuous across the top level
        a1, _, c1 = A5Grid.face_apk(t, s, top - 1e-9); a2, _, c2 = A5Grid.face_apk(t, s, top + 1e-9)
        @test a1 ≈ a2 rtol = 1e-6
        @test c1 ≈ c2 rtol = 1e-6
        @test A5Grid.face_apk(t, s, z0[s]) == (0.0, 0.0, 0.0)
        @test A5Grid.face_apk(t, s, z0[s] - 1.0) == (0.0, 0.0, 0.0)
    end
    for wse in (-1.0, 0.5, 100.0)
        @test A5Grid.face_apk(t, 6, wse) == (0.0, 0.0, 0.0)
    end
end

const _AOI = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
  "coordinates":[[[172.610,-43.515],[172.630,-43.515],[172.630,-43.500],[172.610,-43.500],[172.610,-43.515]]]}}"""
const _A5 = A5Grid.mesh_for_aoi(_AOI, 16)
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

@testset "P2 — per-vertex diamond values match the per-edge lookups" begin
    st = _quiet(() -> initialise_flow_model(_A5, StandardFlow(); gradient_correction = true,
                                            face_flux_method = :diamond))
    T = st.diamond_table; E = st.edges; n = length(st.cell_ids)
    wse = 10 .* rand(n); depth = rand(n) .- 0.3
    veta = diamond_vertex_eta!(Float64[], T, wse)
    vwet = diamond_vertex_wet!(Bool[], T, depth, HFLOW_THRESHOLD)
    @test length(veta) == length(T.vertices) == length(T.vptr) - 1
    @test all(T.vcell[T.vptr[v]:T.vptr[v + 1] - 1] == T.vertices[v].cell_indices &&
              T.vw[T.vptr[v]:T.vptr[v + 1] - 1] == T.vertices[v].weights for v in eachindex(T.vertices))
    nvalid = 0
    for e in 1:E.n_edges
        a = diamond_dWSE_n(T, e, wse[E.cell_i[e]], wse[E.cell_j[e]], wse, E.L[e])
        b = diamond_dWSE_n_v(T, e, wse[E.cell_i[e]], wse[E.cell_j[e]], veta, E.L[e])
        @test a === b || (isnan(a) && isnan(b))
        @test diamond_stencil_wet(T, e, depth, HFLOW_THRESHOLD) ==
              diamond_stencil_wet_v(T, e, vwet)
        nvalid += T.valid[e]
    end
    @test nvalid > 0.9 * E.n_edges
    for v in eachindex(T.vertices)
        x = eta_at_vertex(T.vertices[v], wse)
        @test x === veta[v] || (isnan(x) && isnan(veta[v]))
    end
end

@testset "P3 — velocity on demand; workspace reuse" begin
    for kw in ((; ), (; gradient_correction = true, face_flux_method = :diamond_wetdry,
                        friction_model = :vector, q_centre_mode = :consistent))
        mk() = begin
            s = _quiet(() -> initialise_flow_model(_A5, StandardFlow(); kw...))
            for k in eachindex(s.ghost_cell_bc); s.ghost_cell_bc[k] = Closed; end
            x = [c.center_lon for c in _A5.cells]
            s.elevation .= 1e-3 .* 111_000 .* (x .- minimum(x))
            s.water_depth .= max.(0.0, 0.5 .- s.elevation); s.volume .= s.water_depth .* s.cell_area
            s
        end
        a = mk(); b = mk()
        for _ in 1:20
            step_standard!(a, 1.0)
            step_standard!(b, 1.0; velocity = false)
        end
        _compute_velocity!(b)
        @test a.vel_u == b.vel_u && a.vel_v == b.vel_v && a.velocity == b.velocity
        @test a.edges.flux == b.edges.flux && a.volume == b.volume
        @test length(a.ws.edge_vol) == a.edges.n_edges && length(a.ws.dV) == length(a.cell_ids)
        bytes = @allocated step_standard!(a, 1.0)
        @test bytes < 200_000                              # was ~3 MB/step before WP4
    end
end

println("\nWP4 equivalence tests complete.")
