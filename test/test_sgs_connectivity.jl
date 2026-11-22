# test/test_sgs_connectivity.jl — SGS connectivity and flow-path geometry (0.7.0)
#
#   K1  _cell_connectivity: a wall across a cell's working grid cuts the far
#       side off up to its crest, or to the sill of a gap through it; a wall
#       one point thick and diagonal to the grid holds (4-connectivity); a pit
#       that reaches no side is not the main basin
#   K2  wet_centroid: interpolation of the wetted-area centroid curve
#   K3  _sgs_path: fully wet cells keep the kernel's geometry; cells less than
#       half wet use the line between their wetted centroids; limits
#   K4  build_sgs_tables! on a hex mesh with a ridge narrower than a cell and
#       no gap: with connectivity the faces behind the ridge have raised sills
#       and water held below the crest does not pass; without it, it does
#       (the "hop" of the scheme up to 0.6.0; registry case planar_wall)
#
# Run:  julia --project=. test/test_sgs_connectivity.jl
ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using Test, Logging

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

"A full nx × ny working grid: (z, idx, gi, gj, sidemask) with sides 1 = west, 2 = east, 3 = south, 4 = north."
function _grid(zf, nx, ny)
    idx = zeros(Int32, nx, ny); z = Float64[]; gi = Int32[]; gj = Int32[]; sm = UInt32[]
    for j in 1:ny, i in 1:nx
        push!(z, zf(i, j)); push!(gi, i); push!(gj, j); idx[i, j] = length(z)
        m = UInt32(0)
        i == 1 && (m |= 0x1); i == nx && (m |= 0x2); j == 1 && (m |= 0x4); j == ny && (m |= 0x8)
        push!(sm, m)
    end
    return z, idx, gi, gj, sm
end

@testset "K1 — _cell_connectivity" begin
    nx = ny = 21
    slope(i, j) = 0.01 * i                                   # rises to the east
    # (a) plain slope: one basin, every point connected at its own level
    z, idx, gi, gj, sm = _grid(slope, nx, ny)
    seed, conn = A5Grid._cell_connectivity(z, idx, gi, gj, sm, 2)
    @test z[seed] == minimum(z)
    @test conn == z
    # (b) a wall (crest 3 m) across the middle: the west side is the main basin
    wall(i, j) = i == 11 ? 3.0 : slope(i, j)
    z, idx, gi, gj, sm = _grid(wall, nx, ny)
    seed, conn = A5Grid._cell_connectivity(z, idx, gi, gj, sm, 2)
    @test gi[seed] < 11
    @test all(conn[p] == z[p] for p in eachindex(z) if gi[p] < 11)
    @test all(conn[p] == 3.0 for p in eachindex(z) if gi[p] > 11)   # behind the wall: its crest
    @test all(conn .>= z)
    # (c) a gap through the wall with its sill at 1.5 m
    gap(i, j) = (i == 11 && j == 7) ? 1.5 : wall(i, j)
    z, idx, gi, gj, sm = _grid(gap, nx, ny)
    seed, conn = A5Grid._cell_connectivity(z, idx, gi, gj, sm, 2)
    @test all(conn[p] == 1.5 for p in eachindex(z) if gi[p] > 11)
    @test conn[idx[11, 7]] == 1.5 && conn[idx[11, 8]] == 3.0
    # (d) a wall one point thick, diagonal to the grid: no leak through the corners
    diag(i, j) = i == j ? 3.0 : (i < j ? 0.0 : 0.5)
    z, idx, gi, gj, sm = _grid(diag, nx, ny)
    seed, conn = A5Grid._cell_connectivity(z, idx, gi, gj, sm, 2)
    @test gi[seed] < gj[seed]
    @test all(conn[p] == 3.0 for p in eachindex(z) if gi[p] > gj[p])
    # (e) a pit in the middle, below everything, reaching no side; the plain
    # around it slopes; need = 2. The pit spills at its rim into the plain:
    # the main basin is the first region to reach two sides, and no point of
    # the plain is held above its own level by the pit
    pit(i, j) = (9 <= i <= 12 && 9 <= j <= 12) ? -2.0 : (8 <= i <= 13 && 8 <= j <= 13 ? 1.0 : slope(i, j))
    z, idx, gi, gj, sm = _grid(pit, nx, ny)
    seed, conn = A5Grid._cell_connectivity(z, idx, gi, gj, sm, 2)
    @test count_ones(sm[seed]) >= 1 || z[seed] <= 0.02     # the seed is on the plain's low edge …
    @test z[seed] > -2.0                                      # … not in the pit
    @test all(conn[p] == z[p] for p in eachindex(z) if z[p] >= 0.0 && z[p] < 1.0)   # the plain
    @test all(conn[p] == 1.0 for p in eachindex(z) if z[p] == -2.0)                  # the pit: from its rim
    # (f) empty grid
    @test A5Grid._cell_connectivity(Float64[], zeros(Int32, 0, 0), Int32[], Int32[], UInt32[], 2) == (0, Float64[])
end

@testset "K2 — wet_centroid" begin
    eb = [0.0, 1.0, 2.0]
    t = A5Grid.SGSTable(eb, [0.0, 1.0, 3.0], [0.0, 1.0, 2.0], 2.0, 0.0, 2.0,
                        fill(NaN, 1), zeros(1), zeros(A5Grid.SGS_FACE_LEVELS, 1),
                        zeros(A5Grid.SGS_FACE_LEVELS, 1), zeros(A5Grid.SGS_FACE_LEVELS, 1),
                        [-8.0, -4.0, 0.0], [2.0, 1.0, 0.0])
    @test wet_centroid(t, -1.0) == (-8.0, 2.0)               # dry: the lowest point of the main basin
    @test wet_centroid(t, 0.0) == (-8.0, 2.0)
    @test all(isapprox.(wet_centroid(t, 0.5), (-6.0, 1.5); atol = 1e-9))
    @test all(isapprox.(wet_centroid(t, 1.5), (-2.0, 0.5); atol = 1e-9))
    @test wet_centroid(t, 5.0) == (0.0, 0.0)                 # submerged: the cell's centroid
    # the 11-argument constructor (no centroid curve): the cell centre at every level
    t0 = A5Grid.SGSTable(eb, [0.0, 1.0, 3.0], [0.0, 1.0, 2.0], 2.0, 0.0, 2.0,
                         fill(NaN, 1), zeros(1), zeros(2, 1), zeros(2, 1), zeros(2, 1))
    @test wet_centroid(t0, 0.7) == (0.0, 0.0)
end

@testset "K3 — _sgs_path" begin
    L = 20.0
    # fully wet: the kernel's own geometry, whatever the centroids
    @test _sgs_path(20.0, 0.0, 1.0, 0.0, L, 0.9, 3.0, 4.0, -2.0, 1.0, 1.0, 1.0) == (1.0, L * 0.9)
    # channels (wetted fraction ≤ ½ on average): area projected on, slope along, the centroid line
    a, Lg = _sgs_path(20.0, 0.0, 1.0, 0.0, L, 0.9, 0.0, 0.0, 0.0, 0.0, 0.25, 0.25)
    @test a == 1.0 && Lg == 20.0                              # aligned: nothing changes but cos θ
    # a face cut at 45° by a channel whose water lies 12 m apart
    u = (12 / sqrt(2), 12 / sqrt(2))
    a, Lg = _sgs_path(20.0, 0.0, 1.0, 0.0, L, 1.0, 0.0, 0.0, u[1] - 20.0, u[2], 0.2, 0.3)
    @test a ≈ 1 / sqrt(2) && Lg ≈ 12.0
    # half-way between (mean wetted fraction ¾): the mean of the two
    a2, Lg2 = _sgs_path(20.0, 0.0, 1.0, 0.0, L, 1.0, 0.0, 0.0, u[1] - 20.0, u[2], 0.75, 0.75)
    @test a2 ≈ 0.5 * (1 / sqrt(2) + 1) && Lg2 ≈ 0.5 * (12.0 + 20.0)
    # limits: a ≥ 0.1 (flow along the face), path ≥ 0.2 L (centroids nearly coincident)
    a, Lg = _sgs_path(20.0, 0.0, 1.0, 0.0, L, 1.0, 0.0, 0.0, -20.0, 5.0, 0.1, 0.1)
    @test a == 0.1 && Lg == 5.0
    a, Lg = _sgs_path(20.0, 0.0, 1.0, 0.0, L, 1.0, 0.0, 0.0, -19.0, 0.0, 0.1, 0.1)
    @test a == 1.0 && Lg == 0.2 * L
    # no geometry (NaN centre vector): the kernel's own
    @test _sgs_path(NaN, NaN, 1.0, 0.0, L, 0.8, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1) == (1.0, L * 0.8)
end

# ── K4: a ridge narrower than a cell, no gap ───────────────────────────────
const DIR = mktempdir()
A5Grid._find_python()
const HEX = joinpath(DIR, "hex.parquet")
run(pipeline(Cmd([A5Grid.PYTHON_CMD[], joinpath(@__DIR__, "..", "mesh", "mesh_families.py"),
                  "hex", "--dx", "30", "--bbox", "0.0,0.0,0.003,0.003", "--out", HEX]); stdout = devnull))
# DEM (EPSG:4326, 2 m pixels, 400 × 400 m): a plane rising 0.2% to the east
# with a 2 m-high ridge 6 m wide at x = 176–182 m, from edge to edge.
const DEM = joinpath(DIR, "dem.tif")
const X_RIDGE = (176.0, 182.0); const CREST = 2.0
let px = rad2deg(2.0 / A5Grid._EARTH_R), nx = 200, ny = 200
    z = Float32[0.002 * (2.0 * i) + (X_RIDGE[1] <= 2.0 * i <= X_RIDGE[2] ? CREST : 0.0) for i in 0:(nx-1), j in 0:(ny-1)]
    ArchGDAL.create(DEM; driver = ArchGDAL.getdriver("GTiff"), width = nx, height = ny,
                    nbands = 1, dtype = Float32) do ds
        ArchGDAL.setgeotransform!(ds, [-0.0001, px, 0.0, 0.0031, 0.0, -px])
        ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
        ArchGDAL.write!(ds, z, 1)
    end
end
_mesh(conn) = _quiet() do
    m = load_mesh_geoparquet(HEX)
    sample_dem!(m, DEM)
    build_sgs_tables!(m, DEM; n_bins = 40, connectivity = conn)
    m
end
"x (m, DEM frame) of a longitude"
_x(lon) = (lon + 0.0001) * A5Grid._EARTH_R * π / 180

"Fill the cells west of the ridge to level η, run closed for `nsteps`, return (volume east of the ridge, total)."
function _hop(mesh, η; nsteps = 600, dt = 0.5)
    st = _quiet(() -> initialise_flow_model(mesh, SGSFlow(); friction_model = :vector, limiter = :positivity,
                                            cfl_mode = :local, q_centre_theta = 1.0))
    for k in eachindex(st.ghost_cell_bc); st.ghost_cell_bc[k] = Closed; end
    xs = [_x(c.center_lon) for c in mesh.cells]
    for i in eachindex(st.volume)
        t = st.sgs_tables[i]
        (isnan(t.z_min) || xs[i] > X_RIDGE[1] - 20.0 || t.z_max > η) && continue
        st.volume[i] = t.vol_curve[end] + t.cell_area * (η - t.z_max)
    end
    v0 = sum(st.volume)
    for _ in 1:nsteps; step_sgs!(st, dt); end
    @assert isapprox(sum(st.volume), v0; rtol = 1e-10)
    east = sum(st.volume[i] for i in eachindex(xs) if xs[i] > X_RIDGE[2] + 20.0; init = 0.0)
    return east, v0
end

@testset "K4 — a ridge inside the cells holds water below its crest (connectivity)" begin
    on = _mesh(true); off = _mesh(false)
    @test on.static_vars["sgs_connectivity"][1] == 1.0 && off.static_vars["sgs_connectivity"][1] == 0.0
    @test maximum(off.array_vars["sgs_face_raise"]) == 0.0
    r = on.array_vars["sgs_face_raise"]
    @test count(>(1.0), r) >= 10                             # faces behind the ridge: sill ≈ its crest
    @test maximum(r) <= CREST + 0.5
    # every raised face belongs to a cell on the ridge or to a neighbour of one
    # (a face is shared: the cell east of a ridge cell sees the same raised sill)
    for (i, c) in enumerate(on.cells), k in axes(r, 1)
        r[k, i] > 0.2 && @test X_RIDGE[1] - 65.0 < _x(c.center_lon) < X_RIDGE[2] + 65.0
    end
    # the tables of the two meshes are otherwise the same (storage is level-pool in both)
    @test on.array_vars["sgs_vol_curve"] == off.array_vars["sgs_vol_curve"]
    # water held 0.8 m below the crest (ground at the ridge ≈ 0.35 m, crest ≈ 2.35 m)
    η = 1.5
    east_on, v0 = _hop(on, η)
    east_off, _ = _hop(off, η)
    @test v0 > 1000.0
    @test east_on <= 1e-9 * v0                               # nothing passes
    @test east_off > 0.01 * v0                               # the hop of the scheme up to 0.6.0
    println("  K4: of $(round(Int, v0)) m³ held 0.8 m below the crest, east of the ridge after 300 s: ",
            "$(round(east_on; digits = 3)) m³ with connectivity, $(round(Int, east_off)) m³ without")
end

println("\nSGS connectivity tests complete.")
