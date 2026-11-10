# Compare: run vs reference (compare_fields, compare) and twin-run divergence.

"A format-1 copy of the fixture's frames (FIX_H5 cell order), depths and speeds mapped by `f`."
function _write_v1_frames(path, f::Function)
    r = FA.load_run(FIX_H5)
    h5open(path, "w") do o
        g = create_group(o, "mesh")
        g["cell_ids"] = r.mesh.ids; g["center_lons"] = r.mesh.lon
        g["center_lats"] = r.mesh.lat; g["elevations"] = r.mesh.elevation
        fr = create_group(o, "frames")
        for k in 1:FA.nframes(r)
            gk = create_group(fr, lpad(string(k), 6, '0'))
            gk["t"] = r.t[k]
            gk["water_depth"] = f(FA.depth(r, k)); gk["volume"] = f(FA.volume(r, k))
            for v in ("velocity", "vel_u", "vel_v")
                gk[v] = f(FA.field(r, v, k))
            end
        end
    end
    return FA.load_run(path; parquet = FIX_MESH)
end

@testset "A-compare — compare_fields on known arrays" begin
    h_ref = [0.0, 0.05, 0.2, 0.5, 4.0, NaN]
    h     = [0.02, 0.0, 0.3, 0.0, 4.0, 1.0]
    area  = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    s_ref = [0.0, 0.0, 0.0, 0.0, 2.0, 0.0]; s = [9.0, 9.0, 9.0, 9.0, 1.0, 9.0]
    c = FA.compare_fields(h_ref, h, area; s_ref, s, threshold = 0.1)
    @test c.n_wet == 5                                   # cell 6: NaN in the reference
    @test c.rms ≈ sqrt((0.02^2 + 0.05^2 + 0.1^2 + 0.5^2) / 5)
    @test (c.hits, c.misses, c.false_alarms) == (2, 1, 1)
    @test c.csi == 0.5 && c.bias == 1.0
    @test c.storage_ref ≈ 22.7 && c.storage ≈ 26.92 && c.storage_ratio ≈ 26.92 / 22.7
    @test c.n_deep == 1 && c.speed_ref == 2.0 && c.speed == 1.0 && c.speed_ratio == 0.5
    @test isnan(FA.compare_fields(h_ref, h, area).speed_ratio)          # no speeds given
    @test isnan(FA.compare_fields(zeros(3), zeros(3), ones(3); s_ref = ones(3), s = ones(3)).speed_ratio)
    @test_throws ErrorException FA.compare_fields(h_ref, h[1:5], area)
end

@testset "A-compare — compare(ref, run)" begin
    r = FA.load_run(FIX_H5)
    th = 0.03                                    # the fixture is shallow (≈ 7 cm at most)
    c = FA.compare(r, r; threshold = th)
    @test c.k_ref == c.k_run == FA.nframes(r) && c.t == c.t_ref == 600.0
    @test c.rms == 0.0 && c.csi == 1.0 && c.bias == 1.0 && c.storage_ratio == 1.0
    @test c.n_wet > 0 && c.hits > 0 && c.misses == 0 && c.false_alarms == 0
    @test c.storage == FA.storage(r)
    @test c.n_deep == 0 && isnan(c.speed_ratio)  # no cell ≥ 3 m deep
    @test FA.compare(r, r).hits == 0             # nothing above 0.1 m: CSI and bias 0
    k = FA.frame_at(r, 300.0)
    @test FA.compare(r, r; k_run = k).k_ref == k
    @test FA.compare(r, r; k_ref = k).k_run == k
    # a perturbed copy: depths (and speeds, volumes) × 1.1
    p = _write_v1_frames(joinpath(FIX_DIR, "x1.1.h5"), x -> 1.1 .* x)
    @test p.format == 1
    c = FA.compare(r, p; threshold = th)
    h = FA.depth(r)
    w = filter(i -> isfinite(h[i]) && h[i] > 0.01 / 1.1, eachindex(h))
    @test c.n_wet == length(w)
    @test c.rms ≈ 0.1 * sqrt(sum(abs2, h[w]) / length(w))
    @test c.storage_ratio ≈ 1.1
    @test c.hits == count(>(th), h) > 0 && c.misses == 0
    @test c.false_alarms == count(x -> th / 1.1 < x <= th, h) > 0
    @test c.csi ≈ c.hits / (c.hits + c.false_alarms)
    @test c.bias ≈ count(>(th / 1.1), h) / count(>(th), h)
    @test FA.compare(p, r; threshold = th).storage_ratio ≈ 1 / 1.1
end

@testset "A-compare — divergence and growth_time" begin
    r = FA.load_run(FIX_H5)
    rows = FA.divergence(r, r)
    @test length(rows) == FA.nframes(r) && [x.t for x in rows] == r.t
    @test all(x -> x.rms == 0 && x.mx == 0 && x.n1cm == 0 && x.urms == 0 && x.dV == 0, rows)
    @test rows[end].nwet == count(>(0.05), FA.depth(r))
    @test isnan(FA.growth_time(rows).tau) && FA.growth_time(rows).n == 0
    # a twin with depths + 2 cm in every cell: rms = max = 2 cm in the wet cells
    p = _write_v1_frames(joinpath(FIX_DIR, "plus2cm.h5"), x -> x .+ 0.02)
    d = FA.divergence(r, p)
    @test d[end].rms ≈ 0.02 && d[end].mx ≈ 0.02 && d[end].n1cm == d[end].nwet
    @test d[end].nwet == count(>(0.03), FA.depth(r))
    @test d[end].urms ≈ 0.02
    m = r.mesh                                            # another cell order: refused
    rv = FA.Run(r.path, "rev", r.format, r.t, FA.MeshGeom(reverse(m.ids), m.lon, m.lat, m.elevation,
                m.area, m.polygons, m.neighbours, m.static, m.id_scheme), r.info, r.flow_model)
    @test_throws ErrorException FA.divergence(r, rv)
    # growth_time: exact on a synthetic exponential; only rows with lo < rms < hi count
    syn = [(t = t, rms = 1e-8 * exp(t / 120)) for t in 0.0:60.0:2400.0]
    g = FA.growth_time(syn)
    @test g.tau ≈ 120.0 && g.n == count(x -> 1e-9 < x.rms < 1e-3, syn)
    @test g.n < length(syn)                               # the late rows are excluded
    @test FA.growth_time(syn; hi = 1e-8 * exp(0.6)).n == 2
end
