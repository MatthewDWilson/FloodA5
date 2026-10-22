# test/test_wp4_stage2.jl
# -----------------------
# WP4 stage 2 (2026-10): the restructured step functions against the frozen
# reference copies (surfacewater/reference_steps.jl), and the spatial cell
# order.
#
#   R1  step_standard!/step_sgs! == step_*_reference!, bit for bit, over many
#       steps and option sets (both limiters; SGS on synthetic tables)
#   R2  _positivity_scale_gather! == _positivity_scale!
#   R3  spatial_cell_order is a permutation that brings neighbours together;
#       a run on the permuted mesh agrees with the file-order run to round-off
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_wp4_stage2.jl

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging, Random

Random.seed!(20261023)

const _AOI = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon",
  "coordinates":[[[172.610,-43.515],[172.630,-43.515],[172.630,-43.500],[172.610,-43.500],[172.610,-43.515]]]}}"""
const _A5 = A5Grid.mesh_for_aoi(_AOI, 16)
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

# A wet, sloping, slightly rough state with mixed open/closed boundaries; for
# SGS, synthetic prismatic tables with simple increasing edge curves.
function _wet_state(mesh, method; kw...)
    s = _quiet(() -> initialise_flow_model(mesh, StandardFlow(); kw...))
    for k in eachindex(s.ghost_cell_bc)
        s.ghost_cell_bc[k] = k % 3 == 0 ? Closed : ZeroGradient
    end
    x = [c.center_lon for c in mesh.cells]; y = [c.center_lat for c in mesh.cells]
    s.elevation .= 1e-3 .* 111_000 .* (x .- minimum(x)) .+ 0.05 .* sin.(7e4 .* y)
    s.water_depth .= max.(0.0, 0.6 .- s.elevation) .+ 0.02 .* (eachindex(x) .% 5 .== 0)
    s.volume .= s.water_depth .* s.cell_area
    if method isa SGSFlow
        nb = 8; nsl = s.max_neighbours
        tbls = SGSTable[]
        for i in eachindex(s.cell_ids)
            z = s.elevation[i]; A = s.cell_area[i]
            eb = collect(range(z, z + 2.0; length = nb))
            vc = A .* (eb .- z) .* (1 .+ 0.1 .* (eb .- z))
            ea = repeat(10.0 .* (eb .- z) .+ 0.5, 1, nsl)
            ep = repeat(10.0 .+ 2.0 .* (eb .- z), 1, nsl)
            push!(tbls, SGSTable(eb, vc, fill(A, nb), A, z, z + 2.0, ea, ep))
        end
        s.sgs_tables = tbls
    end
    return s
end

@testset "R1 — restructured steps == frozen reference ($lab)" for (lab, method, kw) in [
        ("std recommended",    StandardFlow(), scheme_options(:recommended, :standard)),
        ("std legacy",         StandardFlow(), scheme_options(:legacy, :standard)),
        ("std positivity",     StandardFlow(), merge(scheme_options(:recommended, :standard), (; limiter = :positivity))),
        ("std legacy limiter", StandardFlow(), merge(scheme_options(:recommended, :standard), (; limiter = :legacy))),
        ("sgs recommended",    SGSFlow(),      scheme_options(:recommended, :sgs)),
        ("sgs positivity",     SGSFlow(),      merge(scheme_options(:recommended, :sgs), (; limiter = :positivity))),
        ("sgs consistent QC",  SGSFlow(),      merge(scheme_options(:recommended, :sgs),
                                                     (; q_centre_theta = 0.9, q_centre_mode = :consistent))),
    ]
    a = _wet_state(_A5, method; kw...); b = _wet_state(_A5, method; kw...)
    sgs = method isa SGSFlow
    for k in 1:60
        dt = 0.8 + 0.05 * (k % 5)
        sgs ? step_sgs!(a, dt) : step_standard!(a, dt)
        sgs ? step_sgs_reference!(b, dt) : step_standard_reference!(b, dt)
    end
    @test sum(abs, a.edges.flux) + sum(abs, a.edges.flux_Q) > 0      # water moved
    @test a.volume == b.volume && a.water_depth == b.water_depth
    @test a.edges.flux == b.edges.flux && a.edges.flux_Q == b.edges.flux_Q
    @test a.vel_u == b.vel_u && a.qvec_u == b.qvec_u && a.vol_removed == b.vol_removed
end

@testset "R2 — threaded positivity gather == serial" begin
    st = _wet_state(_A5, StandardFlow())
    E = st.edges; ne = E.n_edges
    ev = randn(ne) .* 5.0; ev[1:7:end] .= 0.0
    q1 = randn(ne); q2 = copy(q1); ev1 = copy(ev); ev2 = copy(ev)
    vol = abs.(randn(length(st.cell_ids))) .* 3.0
    _positivity_scale!(ev1, q1, vol, E.cell_i, E.cell_j)
    _positivity_scale_gather!(ev2, q2, vol, E, st, _workspace!(st))
    @test ev1 == ev2 && q1 == q2
end

@testset "R3 — spatial cell order" begin
    p = spatial_cell_order(_A5)
    @test sort(p) == collect(eachindex(_A5.cells))
    m2 = permute_mesh(_A5, p)
    @test [c.id for c in m2.cells] == [c.id for c in _A5.cells][p]
    # initial state from coordinates only (not cell or ghost-edge indices)
    function ordfree(mesh)
        # recommended options: consistent Q-centring is orientation-independent
        # (legacy Q-centring is not, which is why run_flood_model keeps the
        # file order for it)
        s = _quiet(() -> initialise_flow_model(mesh, StandardFlow(); scheme_options(:recommended, :standard)...))
        fill!(s.ghost_cell_bc, ZeroGradient)
        x = [c.center_lon for c in mesh.cells]; y = [c.center_lat for c in mesh.cells]
        s.elevation .= 1e-3 .* 111_000 .* (x .- 172.61) .+ 0.05 .* sin.(7e4 .* y)
        s.water_depth .= max.(0.0, 0.6 .- s.elevation)
        s.volume .= s.water_depth .* s.cell_area
        s
    end
    a = ordfree(_A5); b = ordfree(m2)
    spread(st) = sum(abs.(st.edges.cell_i[1:st.edges.n_edges] .- st.edges.cell_j[1:st.edges.n_edges])) /
                 st.edges.n_edges
    @test spread(b) < 0.5 * spread(a)                                 # neighbours closer in memory
    for _ in 1:30; step_standard!(a, 1.0); step_standard!(b, 1.0); end
    @test maximum(abs, a.volume[p] .- b.volume) <= 1e-9 * maximum(a.volume)
    @test sum(a.volume) ≈ sum(b.volume) rtol = 1e-12
end

# ---------------------------------------------------------------------------
# R4 (stage 3) — active-edge list mode == full pass, on a partly wet domain
# where cells wet, drain and dry again. Mass balance is part of the check:
# total stored + removed volume must match the reference bit for bit.
# (Mutation check, 2026-10-23: dropping the "visit just-dried edges once more"
# rule makes the legacy-limiter configurations fail here. With the hybrid or
# positivity limiter the positivity scaling, which runs over every edge,
# zeroes a stale transfer from an empty donor — so those are self-correcting.)
function _blob_state(method; consistent_sill::Bool, kw...)
    mesh = deepcopy(_A5)
    x = [c.center_lon for c in mesh.cells]; y = [c.center_lat for c in mesh.cells]
    z = 2e-3 .* 111_000 .* (maximum(x) .- x) .+ 0.05 .* sin.(7e4 .* y)   # slopes towards max lon
    consistent_sill && (mesh.static_vars["elevation"] = z)               # sills = max(bed) at init
    s = _quiet(() -> initialise_flow_model(mesh, StandardFlow(); kw...))
    fill!(s.ghost_cell_bc, ZeroGradient)
    s.elevation .= z
    top = sortperm(x)[1:12]                                             # wet blob, high side
    fill!(s.volume, 0.0); fill!(s.water_depth, 0.0)
    s.water_depth[top] .= 0.4; s.volume[top] .= 0.4 .* s.cell_area[top]
    if method isa SGSFlow
        nb = 8; nsl = s.max_neighbours
        s.sgs_tables = [begin
            zi = z[i]; A = s.cell_area[i]; eb = collect(range(zi, zi + 2.0; length = nb))
            SGSTable(eb, A .* (eb .- zi) .* (1 .+ 0.1 .* (eb .- zi)), fill(A, nb), A, zi, zi + 2.0,
                     repeat(10.0 .* (eb .- zi) .+ 0.5, 1, nsl), repeat(10.0 .+ 2.0 .* (eb .- zi), 1, nsl))
        end for i in eachindex(s.cell_ids)]
    end
    return s, top[1]
end

@testset "R4 — active-edge list == full pass ($lab, sills $(cs ? "consistent" : "inconsistent"))" for
        (lab, method, kw) in [
            ("std recommended", StandardFlow(), scheme_options(:recommended, :standard)),
            ("std legacy",      StandardFlow(), scheme_options(:legacy, :standard)),
            ("std positivity",  StandardFlow(), merge(scheme_options(:recommended, :standard), (; limiter = :positivity))),
            ("sgs recommended", SGSFlow(),      scheme_options(:recommended, :sgs)),
            ("sgs positivity",  SGSFlow(),      merge(scheme_options(:recommended, :sgs), (; limiter = :positivity))),
        ], cs in (true, false)
    # (a test driver may restrict to one configuration)
    isdefined(Main, :ONLY_CFG) && (lab != Main.ONLY_CFG || cs != Main.ONLY_CS) && continue
    a, src = _blob_state(method; consistent_sill = cs, kw...)
    ACTIVE_MAX_FRACTION[] = 1.0             # list mode whenever possible
    targeted_in_list = false
    b, _   = _blob_state(method; consistent_sill = cs, kw...)
    sgs = method isa SGSFlow
    n_list = 0; min_act = typemax(Int); dried = 0
    for k in 1:400
        dt = 0.8 + 0.05 * (k % 5)
        k <= 250 && (a.volume[src] += 0.2 * dt; b.volume[src] += 0.2 * dt)   # point inflow
        if k in (150, 300)          # drying event: shallow (< 5 mm, i.e. flowing) cells empty
            for s in (a, b), i in eachindex(s.volume)
                if 0.0 < s.volume[i] < 5e-3 * s.cell_area[i]
                    s.volume[i] = 0.0
                    s === a && (dried += 1)           # wet cells emptied by the event
                end
            end
        end
        if k == 200                 # targeted: empty both cells of flowing edges
            fl = sgs ? a.edges.flux_Q : a.edges.flux
            for e in findall(!=(0.0), fl[1:a.edges.n_edges])[1:min(5, end)]
                for s in (a, b), c in (a.edges.cell_i[e], a.edges.cell_j[e])
                    s.volume[c] = 0.0
                end
                dried += 2
            end
        end
        if !sgs                                     # run_simulation!'s depth sync
            for s in (a, b); s.water_depth .= s.volume ./ s.cell_area; end
        end
        sgs ? step_sgs!(a, dt) : step_standard!(a, dt)
        sgs ? step_sgs_reference!(b, dt) : step_standard_reference!(b, dt)
        k == 201 && (targeted_in_list = a.ws.n_act > 0)
        if a.ws.n_act > 0 && k > 2
            n_list += 1; min_act = min(min_act, a.ws.n_act)
        end
    end
    @test a.volume == b.volume && a.water_depth == b.water_depth
    @test a.edges.flux == b.edges.flux && a.edges.flux_Q == b.edges.flux_Q
    @test a.qvec_u == b.qvec_u && a.vol_removed == b.vol_removed
    @test sum(a.volume) + a.vol_removed == sum(b.volume) + b.vol_removed     # mass balance
    @test n_list > 100                                      # list mode was exercised
    # …on a genuinely partial set when sills = max(bed) (inconsistent sills make
    # every edge "dry-live", which is safe but saves nothing)
    cs && @test min_act < a.edges.n_edges ÷ 2
    @test a.vol_removed > 0                                 # water reached the boundary
    @test dried > 0                                         # wet cells emptied → just-dried path
    cs && @test targeted_in_list                            # …and handled in list mode
    ACTIVE_MAX_FRACTION[] = 0.5
end

println("\nWP4 stage 2 tests complete.")
