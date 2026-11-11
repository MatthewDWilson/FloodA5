# Runs: loading, frames, fields (format 2 fixture).
@testset "A-runs — load_run on a format-2 file" begin
    r = FA.load_run(FIX_H5)
    @test r.format == 2 && r.flow_model == "standard"
    @test FA.nframes(r) == FO.n_frames(FIX_H5) && FA.times(r)[end] == 600.0
    n = length(r.mesh)
    @test length(r.mesh.area) == n && length(r.mesh.polygons) == n && size(r.mesh.neighbours, 2) == n
    @test FA.frame_at(r, r.t[4] - 5.0) == 4
    @test_throws ErrorException FA.frame_at(r, 10_000.0; tol = 60)
    @test FA.matched_frames(r, r) == [(k, k) for k in 1:FA.nframes(r)]
    h = FA.depth(r)
    @test isequal(FA.wse(r), r.mesh.elevation .+ h)
    @test FA.storage(r) ≈ FA.stored_volume(r) rtol = 1e-9
    @test FA.storage(r) ≈ 2.0 * 600.0 rtol = 0.01               # closed domain, 2 m³/s × 600 s
    i = FA.source_cell(r.mesh, FIX_SRC...)
    x, y = FA.local_xy(r.mesh, r.mesh.lon[i], r.mesh.lat[i])
    V1 = replace(FA.volume(r, FA.frame_at(r, 60.0)), NaN => 0.0)
    @test hypot(sum(V1 .* x), sum(V1 .* y)) / sum(V1) < 12.5   # early volume centroid in the source cell
    @test x[i] == 0.0 && y[i] == 0.0
end

@testset "A-runs — storage of SGS runs is Σ volume" begin
    r = FA.load_run(FIX_H5)
    s = FA.Run(r.path, r.label, r.format, r.t, r.mesh, r.info, "sgs")   # the same file, read as SGS
    @test FA.storage(s) == FA.stored_volume(s)
    @test FA.compare(s, s).storage_ratio == 1.0 && FA.compare(r, s).storage == FA.stored_volume(s)
end
