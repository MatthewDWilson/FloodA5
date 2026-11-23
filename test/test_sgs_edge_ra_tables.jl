"""
test_sgs_edge_ra_tables.jl — unit tests of the SGS face tables (layout 3, 0.7.0).

One table per face, built from the DEM profile along the face
(`A5Grid.face_table_from_profile`) and read by `A5Grid.face_apk`: flow area A,
wetted width P and conveyance integral K = Σ dx·h^(5/3) on the face's own
levels, continued analytically above the tabulated range, and (layout 4,
0.8.0) the critical integral C = Σ dx·h^(3/2) read by `A5Grid.face_crit` for
the critical-flow limit (`_sgs_critical_q`). The tests compare
them with analytic cross-sections (no DEM needed), including far above the
highest point of the profile — where the tables before 0.7.0 stopped growing
(case sgs_strip) — and the two conveyance forms of `_sgs_radius`.

(The file keeps its name from the R–A edge curves it tested until 0.7.0.)

Run:
    julia --project=. test/test_sgs_edge_ra_tables.jl
"""

using Test, Printf

if abspath(PROGRAM_FILE) == @__FILE__

    project_dir      = dirname(dirname(abspath(@__FILE__)))
    flood_model_path = joinpath(project_dir, "FloodModel.jl")
    isfile(flood_model_path) || error("FloodModel.jl not found at $flood_model_path")

    ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
    include(flood_model_path)
    delete!(ENV, "FLOODMODEL_INCLUDE_ONLY")

    println("=" ^ 62)
    println("SGS face table unit tests")
    println("=" ^ 62)

    "An SGSTable whose slot `slot` (of 5) holds the face table of `profile` (width W)."
    function _table(profile::Vector{Float64}, W::Float64, slot::Int; act = profile)
        nf = A5Grid.SGS_FACE_LEVELS
        z0 = fill(NaN, 5); zr = zeros(5); A = zeros(nf, 5); P = zeros(nf, 5); K = zeros(nf, 5)
        C = zeros(nf, 5)
        z0[slot], zr[slot], A[:, slot], P[:, slot], K[:, slot], C[:, slot] =
            A5Grid.face_table_from_profile(profile, W / length(profile); act = act)
        return A5Grid.SGSTable([0.0, 1.0], [0.0, 1.0], [0.0, 1.0], 1000.0, 0.0, 1.0, z0, zr, A, P, K,
                               zeros(2), zeros(2), C)
    end
    "Exact critical integral Σ dx·h^(3/2) of a profile at level η."
    _exact_c(profile, W, η) = sum(max(η - z, 0.0)^1.5 for z in profile) * W / length(profile)
    "Exact (A, P, K) of a profile at level η."
    function _exact(profile, W, η)
        dx = W / length(profile)
        (sum(max(η - z, 0.0) for z in profile) * dx, count(<=(η), profile) * dx,
         sum(max(η - z, 0.0)^(5 / 3) for z in profile) * dx)
    end

    @testset "F1 — flat face: A and P exact; K linear between levels (≤ 5% at 1 mm, 0.1% from 0.1 m)" begin
        W = 10.0; z = 5.0
        t = _table(fill(z, 20), W, 1)
        for d in (1e-3, 0.1, 0.5, 1.0, 2.0, 10.0)
            A, P, K = A5Grid.face_apk(t, 1, z + d)
            @test A ≈ d * W rtol = 1e-12
            @test P ≈ W rtol = 1e-12
            @test K ≈ W * d^(5 / 3) rtol = (d < 0.05 ? 0.05 : 1e-3)
        end
    end

    @testset "F2 — dry at and below the sill; unused slots are dry" begin
        t = _table(fill(3.0, 16), 8.0, 2)
        @test A5Grid.face_apk(t, 2, 3.0) == (0.0, 0.0, 0.0)
        @test A5Grid.face_apk(t, 2, 2.0) == (0.0, 0.0, 0.0)
        @test A5Grid.face_apk(t, 2, 3.5)[1] > 0.0
        for s in (1, 3, 4, 5)
            @test A5Grid.face_apk(t, s, 100.0) == (0.0, 0.0, 0.0)
        end
        @test isnan(t.face_z0[1]) && t.face_z0[2] == 3.0
    end

    # A channel a quarter of the face wide and 1 m deep in a flat floodplain:
    # the strip case's cross-channel face (sgs_strip).
    chan = vcat(fill(10.0, 300), fill(9.0, 200), fill(10.0, 300)); Wc = 20.0

    @testset "F3 — channel in a floodplain: A, P, K against the profile" begin
        t = _table(chan, Wc, 3)
        @test t.face_z0[3] == 9.0
        @test t.face_zr[3] ≈ 3.0                       # relief 1 m + max(2·relief, 0.5)
        for η in (9.05, 9.3, 9.73, 9.99, 10.01, 10.37, 11.0, 12.0)   # tabulated range
            A, P, K = A5Grid.face_apk(t, 3, η); Ae, Pe, Ke = _exact(chan, Wc, η)
            # at the bank the width jumps from 5 to 20 m inside one level interval (5 cm there)
            @test A ≈ Ae rtol = (abs(η - 10.0) < 0.05 ? 0.03 : 0.01)
            @test K ≈ Ke rtol = 0.02
            η > 10.02 && @test P ≈ Pe rtol = 1e-9
            η < 9.99 && @test P ≈ Pe rtol = 1e-9
        end
        for η in (12.5, 15.0, 30.0, 100.0)             # above it: the continuation
            A, P, K = A5Grid.face_apk(t, 3, η); Ae, Pe, Ke = _exact(chan, Wc, η)
            @test A ≈ Ae rtol = 1e-9                   # area is exact (full width)
            @test P ≈ Wc rtol = 1e-12
            @test K ≈ Ke rtol = 0.02
        end
    end

    @testset "F4 — nothing stops growing above the highest point" begin
        t = _table(chan, Wc, 3)
        Ap, Kp = 0.0, 0.0
        for η in range(9.0, 40.0; length = 400)
            A, P, K = A5Grid.face_apk(t, 3, η)
            @test A >= Ap && K >= Kp
            Ap, Kp = A, K
        end
        A1, _, K1 = A5Grid.face_apk(t, 3, 14.0); A2, _, K2 = A5Grid.face_apk(t, 3, 15.0)
        @test A2 - A1 ≈ Wc rtol = 1e-9                 # dA/dη = full width, 4 m above the floodplain
        @test K2 > 1.3K1
    end

    @testset "F5 — V section: A within 1%, K within 2% (3% and 6% at 5 cm depth)" begin
        n = 400; W = 12.0
        prof = [2.0 + abs(i - (n + 1) / 2) / (n / 2) * 1.5 for i in 1:n]
        t = _table(prof, W, 1)
        for η in (2.05, 2.5, 3.0, 3.5, 4.0, 8.0, 20.0)
            A, P, K = A5Grid.face_apk(t, 1, η); Ae, Pe, Ke = _exact(prof, W, η)
            sh = η < 2.1                              # 5 cm of water in the V: A = 0.01 m²
            @test A ≈ Ae rtol = (sh ? 0.03 : 0.01)
            @test K ≈ Ke rtol = (sh ? 0.06 : 0.02)
        end
    end

    @testset "F6 — conveyance forms (_sgs_radius)" begin
        t = _table(chan, Wc, 3)
        # in bank the two forms agree (a rectangle: every strip has the same depth)
        A, P, K = A5Grid.face_apk(t, 3, 9.73)
        @test _sgs_radius(A, P, K, true) ≈ _sgs_radius(A, P, K, false) rtol = 0.01
        # overbank: 2-D form A·R^(2/3) = K; lumped A^(5/3)/P^(2/3) is smaller
        A, P, K = A5Grid.face_apk(t, 3, 10.4)
        R2 = _sgs_radius(A, P, K, true); Rl = _sgs_radius(A, P, K, false)
        @test A * R2^(2 / 3) ≈ K rtol = 1e-12
        @test Rl ≈ A / P rtol = 1e-12
        @test A * Rl^(2 / 3) < 0.9K
        @test _sgs_radius(0.0, 0.0, 0.0, true) == 0.0 && _sgs_radius(0.0, 0.0, 0.0, false) == 0.0
    end

    @testset "F7 — activation levels: a point conveys only from its level" begin
        # the floodplain half right of the channel is behind a bank at 10.8 m
        act = vcat(fill(10.0, 300), fill(9.0, 200), fill(10.8, 300))
        t = _table(chan, Wc, 1; act = act)
        A, P, K = A5Grid.face_apk(t, 1, 10.5)
        @test P ≈ 12.5 rtol = 0.02                     # channel + left floodplain only
        A, P, K = A5Grid.face_apk(t, 1, 11.5); Ae, Pe, Ke = _exact(chan, Wc, 11.5)
        @test P ≈ Wc rtol = 1e-9
        @test A ≈ Ae rtol = 0.01
        @test K ≈ Ke rtol = 0.02
    end

    @testset "F8 — fixture constructor (per-slot curves on the cell's levels)" begin
        eb = collect(range(5.0, 8.0; length = 100)); W = 10.0
        ea = zeros(100, 5); ep = zeros(100, 5)
        ea[:, 2] .= W .* (eb .- 5.0); ep[:, 2] .= W; ep[1, 2] = 0.0
        t = A5Grid.SGSTable(eb, zeros(100), zeros(100), 1000.0, 5.0, 8.0, ea, ep)
        for d in (0.2, 1.0, 2.9, 6.0)
            A, P, K = A5Grid.face_apk(t, 2, 5.0 + d)
            @test A ≈ W * d rtol = 0.01
            @test K ≈ W * d^(5 / 3) rtol = 0.02
        end
        @test A5Grid.face_apk(t, 1, 7.0) == (0.0, 0.0, 0.0)
        # no profile: the uniform-depth critical integral P·(A/P)^(3/2)
        @test A5Grid.face_crit(t, 2, 6.0) ≈ W * 1.0^1.5 rtol = 0.02
    end

    @testset "F9 — critical integral and critical discharge (0.8.0)" begin
        # flat face: C = W·H^(3/2); the critical discharge is the
        # broad-crested-weir discharge √g (2/3)^(3/2) W H^(3/2) = 1.705 W H^(3/2)
        W = 4.0; z = 5.0
        t = _table(fill(z, 16), W, 3)
        for H in (0.1, 0.5, 1.07, 3.0, 12.0)
            @test A5Grid.face_crit(t, 3, z + H) ≈ W * H^1.5 rtol = 2e-3
        end
        @test _sgs_critical_q(A5Grid.face_crit(t, 3, z + 1.07)) ≈ 1.7047 * W * 1.07^1.5 rtol = 2e-3
        @test 7.4 < _sgs_critical_q(A5Grid.face_crit(t, 3, z + 1.07)) < 7.6   # the notch of planar_embankment
        @test A5Grid.face_crit(t, 3, z) == 0.0 && A5Grid.face_crit(t, 3, z - 1.0) == 0.0
        @test A5Grid.face_crit(t, 1, z + 1.0) == 0.0                          # unused slot
        # channel in a floodplain: within 2% of the profile at every level
        prof = [i in 9:12 ? 9.0 : 10.0 + 0.01 * i for i in 1:40]; Wc = 20.0
        tc = _table(prof, Wc, 2)
        for η in (9.05, 9.5, 10.0, 10.2, 11.0, 14.0, 30.0)
            @test A5Grid.face_crit(tc, 2, η) ≈ _exact_c(prof, Wc, η) rtol = 0.02
        end
        # activation levels: a point passes critical flow only from its level
        act = copy(prof); act[1:20] .= 12.0
        ta = _table(prof, Wc, 2; act = act)
        @test A5Grid.face_crit(ta, 2, 11.5) ≈ _exact_c(prof[21:end], Wc / 2, 11.5) rtol = 0.02
        # the critical integral grows with the level (the limit never closes a
        # face that is wetter)
        cs = [A5Grid.face_crit(tc, 2, η) for η in range(8.9, 40.0; length = 400)]
        @test all(diff(cs) .>= 0.0)
    end

    println()
    println("=" ^ 62)
    println("All SGS face table unit tests passed.")
    println("=" ^ 62)
end
