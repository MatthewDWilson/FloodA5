# Noise / checkerboard metrics (analysis/noise.jl) on the fixture's square grid.
@testset "A-noise — checkerboard and velocity-noise metrics" begin
    r = FA.load_run(FIX_H5)
    nb = r.mesh.neighbours
    n = length(r.mesh)
    @test size(nb, 1) == 4                                   # Cartesian: 4-connected
    deg = [count(!=(0), view(nb, :, i)) for i in 1:n]
    x, y = FA.local_xy(r.mesh, r.mesh.lon[1], r.mesh.lat[1])
    dx = sqrt(sum(r.mesh.area) / n)                          # 25 m
    @test dx ≈ 25.0 rtol = 1e-3
    lists = [Int[j for j in view(nb, :, i) if j != 0] for i in 1:n]
    zero_v = zeros(n)

    # Smooth paraboloid c(x² + y²): on a square grid r = −c·dx² exactly (4 neighbours).
    # Dry boundary cells, so only cells with 4 wet neighbours count.
    c = 1e-5
    w = 10.0 .+ c .* (x .^ 2 .+ y .^ 2)
    h = [d == 4 ? 1.0 : 0.0 for d in deg]
    m = FA.noise_metrics(w, h, zero_v, zero_v, nb)
    ninner = count(i -> h[i] > 0 && all(j -> h[j] > 0, lists[i]), 1:n)
    @test m.n == ninner > 0
    @test m.r_rms ≈ 1e3 * c * dx^2 rtol = 0.01               # 6.25 mm
    @test m.extrema <= 1 / m.n                               # at most the one minimum
    @test m.dv_rms == 0.0
    @test FA.noise_metrics(w, h, zero_v, zero_v, lists) == m  # list form = matrix form

    # ±a checkerboard (two-colour the grid graph): r = ±2a, every cell a strict extremum.
    col = zeros(Int, n); col[1] = 1; queue = [1]
    while !isempty(queue)
        i = popfirst!(queue)
        for j in lists[i]
            col[j] == 0 && (col[j] = -col[i]; push!(queue, j))
        end
    end
    @test all(i -> all(j -> col[j] == -col[i], lists[i]), 1:n)   # bipartite
    a = 0.01
    hc = ones(n)
    mc = FA.noise_metrics(10.0 .+ a .* col, hc, zero_v, zero_v, nb)
    @test mc.n == count(>=(3), deg)                          # corners (2 nbrs) excluded
    @test mc.r_rms ≈ 2e3 * a rtol = 1e-9                     # 20 mm
    @test mc.extrema == 1.0
    # amplitude below rmin/2: same r_rms, no counted extrema
    @test FA.noise_metrics(10.0 .+ 5e-4 .* col, hc, zero_v, zero_v, nb).extrema == 0.0
    # a checkerboard on the western half only → extrema ≈ 0.5
    half = [col[i] == 1 && x[i] < (minimum(x) + maximum(x)) / 2 for i in 1:n]
    mh = FA.noise_metrics(10.0 .+ a .* half, hc, zero_v, zero_v, nb)
    @test 0.4 < mh.extrema < 0.6

    # Constant velocity → dv_rms 0; a velocity checkerboard → dv_rms = 2|Δv|
    @test FA.noise_metrics(w, hc, fill(0.7, n), fill(-0.2, n), nb).dv_rms ≈ 0 atol = 1e-12
    @test FA.noise_metrics(w, hc, 0.1 .* col, zero_v, nb).dv_rms ≈ 0.2 rtol = 1e-9
    # Nothing wet → NaN
    mz = FA.noise_metrics(w, zeros(n), zero_v, zero_v, nb)
    @test mz.n == 0 && isnan(mz.r_rms) && isnan(mz.extrema) && isnan(mz.dv_rms)
    @test_throws ErrorException FA.noise_metrics(w[1:end-1], hc, zero_v, zero_v, nb)

    # Run level: frame k = the pure function on that frame; series and summary.
    # The fixture has no DEM (elevations NaN, so every WSE is NaN): give it a flat bed.
    @test FA.noise_metrics(r).n == 0
    m0 = r.mesh
    flat = FA.MeshGeom(m0.ids, m0.lon, m0.lat, zeros(n), m0.area, m0.polygons, m0.neighbours,
                       Dict("sgs_z_min" => fill(-0.0, n)), m0.id_scheme)
    rf = FA.Run(r.path, r.label, r.format, r.t, flat, r.info, r.flow_model)
    k = FA.nframes(rf)
    hk = FA.depth(rf, k)
    hmin = 0.005                                             # ≤ 7 cm of water at the end
    uk, vk = FA.field(rf, "vel_u", k), FA.field(rf, "vel_v", k)
    mk = FA.noise_metrics(rf, k; hmin)
    @test mk.n > 0 && isequal(mk, FA.noise_metrics(hk, hk, uk, vk, nb; hmin))
    # SGS: WSE from sgs_z_min (= the bed here: same metrics)
    @test isequal(FA.noise_metrics(rf, k; hmin, sgs = true), mk)
    @test_throws ErrorException FA.noise_metrics(r, k; sgs = true)   # no sgs_z_min
    s = FA.noise_series(rf; every = 4, hmin)
    @test [e.k for e in s] == [kk for kk in 1:k if kk % 4 == 0 || kk == k]
    @test all(e -> e.t == rf.t[e.k], s) && isequal(s[end].m, mk)
    sm = FA.noise_summary(s)
    @test sm.last === s[end]
    valid = filter(e -> e.m.n > 0, s)
    @test !isempty(valid)
    @test sm.max_r_rms == maximum(e.m.r_rms for e in valid)
    @test sm.t_max_r_rms == valid[argmax([e.m.r_rms for e in valid])].t
    @test sm.max_extrema == maximum(e.m.extrema for e in valid)
    @test sm.max_dv_rms == maximum(e.m.dv_rms for e in valid)
    se = FA.noise_summary([(t = 0.0, k = 1, m = mz)])
    @test isnan(se.max_r_rms) && isnan(se.t_max_r_rms) && isnan(se.max_dv_rms)
end
