# test/bench_swe_dambreak.jl
# --------------------------
# R1 (2026-10): polygon consistency of the shallow-water reference solver.
# A circular dam break on a flat, frictionless bed, run on an A5 mesh and on
# a Cartesian grid of the same cell area, and on a fine Cartesian grid (the
# reference profile). Reports, at time T, the radially binned depth profile
# of each run, its difference from the fine profile, and the azimuthal
# spread within each radius band (0 for a perfectly isotropic scheme).
#
#   FLOODA5_NO_VIS=1 julia --project=. -t 8 test/bench_swe_dambreak.jl [--order 2] [--T 60]

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Logging, Printf, Statistics

_arg(name, default) = (i = findfirst(==(name), ARGS); i === nothing ? default : ARGS[i + 1])
const ORDER = parse(Int, _arg("--order", "2"))
const T     = parse(Float64, _arg("--T", "60"))
const R0, H_IN, H_OUT = 400.0, 2.0, 0.5            # dam radius (m), depths inside/outside
const RMAX = 1000.0                                # analysed radius
_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Error))
const _R = A5Grid._EARTH_R

function run_case(mesh)
    mesh.static_vars["elevation"] = zeros(length(mesh.cells))
    s = _quiet(() -> initialise_flow_model(mesh, SWEFlow(); swe_order = ORDER, manning_n = 0.0))
    fill!(s.ghost_cell_bc, Closed)
    lon0 = mean(s.cell_lons); lat0 = mean(s.cell_lats)
    x = deg2rad.(s.cell_lons .- lon0) .* _R .* cosd(lat0); y = deg2rad.(s.cell_lats .- lat0) .* _R
    r = hypot.(x, y)
    s.volume .= [ri < R0 ? H_IN : H_OUT for ri in r] .* s.cell_area
    s.water_depth .= s.volume ./ s.cell_area
    t = 0.0
    while t < T
        dt = min(_cfl_dt(s, SWEFlow(); courant = 0.5), T - t); step_swe!(s, dt; courant = 0.5); t += dt
    end
    return r, s.volume ./ s.cell_area, mean(s.cell_area)
end

function profile(r, h, edges)
    m = fill(NaN, length(edges) - 1); sd = fill(NaN, length(edges) - 1)
    for k in 1:length(edges) - 1
        sel = [h[i] for i in eachindex(r) if edges[k] <= r[i] < edges[k + 1]]
        length(sel) >= 3 || continue
        m[k] = mean(sel); sd[k] = std(sel)
    end
    return m, sd
end

aoi_half = 0.9 * rad2deg(1.6 * RMAX / _R)
aoi = """{"type":"Feature","properties":{},"geometry":{"type":"Polygon","coordinates":[[[$(172.62 - aoi_half / cosd(43.5)),$(-43.5 - aoi_half)],[$(172.62 + aoi_half / cosd(43.5)),$(-43.5 - aoi_half)],[$(172.62 + aoi_half / cosd(43.5)),$(-43.5 + aoi_half)],[$(172.62 - aoi_half / cosd(43.5)),$(-43.5 + aoi_half)],[$(172.62 - aoi_half / cosd(43.5)),$(-43.5 - aoi_half)]]]}}"""
a5 = _quiet(() -> A5Grid.mesh_for_aoi(aoi, 16))
ra5, ha5, area = run_case(a5)
dxm = sqrt(area)                                   # Cartesian grid of the same cell area
cart(dx) = (L = 3.2 * RMAX; d = rad2deg(dx / _R);
            generate_cartesian_mesh_bounds(172.62, -43.5, 172.62 + (L / dx) * d / cosd(-43.5), -43.5 + (L / dx) * d, dx; lat_ref = -43.5))
rc, hc, _ = run_case(cart(dxm))
rf, hf, _ = run_case(cart(dxm / 4))

edges = collect(0.0:dxm:RMAX)
mf, _ = profile(rf, hf, edges)
println(@sprintf("Circular dam break, T = %.0f s, order %d; A5 res-16 (%.0f m² cells, %d cells) vs Cartesian %.1f m and fine %.1f m",
                 T, ORDER, area, length(ha5), dxm, dxm / 4))
@printf("%-10s %10s %10s %10s %10s %10s %10s\n", "r (m)", "fine h", "A5 h", "Cart h", "A5 spread", "Cart spread", "")
ma, sa = profile(ra5, ha5, edges); mc, sc = profile(rc, hc, edges)
for k in eachindex(mf)
    isfinite(mf[k]) || continue
    @printf("%-10.0f %10.4f %10.4f %10.4f %10.4f %10.4f\n", 0.5 * (edges[k] + edges[k + 1]), mf[k], ma[k], mc[k], sa[k], sc[k])
end
ok(v) = filter(isfinite, v)
println(@sprintf("mean |h − fine|: A5 %.4f m, Cartesian %.4f m;  mean azimuthal spread: A5 %.4f m, Cartesian %.4f m",
                 mean(ok(abs.(ma .- mf))), mean(ok(abs.(mc .- mf))), mean(ok(sa)), mean(ok(sc))))
