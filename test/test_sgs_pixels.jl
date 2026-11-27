# test/test_sgs_pixels.jl — SGS tables from the DEM's pixels (0.9.2)
#
# `build_sgs_tables!(mesh, dem; sampling = :pixels)` (`--sgs-sampling pixels`)
# on a Cartesian mesh whose cells are whole k × k blocks of the DEM's pixels
# builds the tables from the pixels themselves, so that an SGS model on the
# coarse mesh stores, at every level, what a model on the DEM's own grid
# stores (HYDRAULICS §15.16).
#
#   X1  cases/scripts/block_dem.jl: k × k block means from the raster's own
#       origin, whole blocks only; the mesh from it lines up with the pixels
#   X2  each cell's volume table is exact: at any level η the stored volume is
#       Σ_p (A/k²)·max(η − z_p, 0) over the block's pixels (and the inverse,
#       wse_from_volume, returns that η); the full area is the cell's area
#   X3  faces: without connectivity a face's sill is the lowest of
#       max(z_in, z_out) over the pixel pairs across it, and both cells hold
#       the same table; a boundary side's sill is its lowest inside pixel
#   X4  a mesh not aligned with the DEM's pixels is refused
#
# Run:  julia --project=. test/test_sgs_pixels.jl
include(joinpath(@__DIR__, "..", "mesh", "A5Grid.jl"))
using .A5Grid
using Test, ArchGDAL, Logging, Random
isdefined(Main, :generate_cartesian_mesh) || include(joinpath(@__DIR__, "..", "mesh", "CartesianGrid.jl"))

_quiet(f) = with_logger(f, SimpleLogger(stderr, Logging.Warn))

const DIR = mktempdir()
const P = 5.0                                  # pixel (m), British National Grid
const X0, Y0 = 338495.0, 557760.0              # north-west corner: not a multiple of 25 m
const NX, NY, K = 32, 21, 5                    # 6 × 4 whole blocks; 2 columns and 1 row dropped
const ZP = let rng = MersenneTwister(7)
    [10.0 + 0.02 * c - 0.03 * r + 1.5 * rand(rng) + (c in 11:12 ? -2.0 : 0.0) for c in 0:(NX-1), r in 0:(NY-1)]
end

function write_tif(path, z, x0, y0, px)
    ArchGDAL.create(path; driver = ArchGDAL.getdriver("GTiff"), width = size(z, 1), height = size(z, 2),
                    nbands = 1, dtype = Float64) do ds
        ArchGDAL.setgeotransform!(ds, [x0, px, 0.0, y0, 0.0, -px])
        ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(27700)))
        ArchGDAL.write!(ds, z, 1)
    end
    return path
end
const FINE = write_tif(joinpath(DIR, "fine.tif"), ZP, X0, Y0, P)
const COARSE = joinpath(DIR, "coarse.tif")

# X1: the block-average script
empty!(ARGS); append!(ARGS, [FINE, COARSE, string(K)])
redirect_stdout(devnull) do
    Base.include(Module(), joinpath(@__DIR__, "..", "cases", "scripts", "block_dem.jl"))
end
empty!(ARGS)

const MESH = _quiet(() -> generate_cartesian_mesh_raster(COARSE))
const NC = length(MESH.cells)
# block (0-based column, row of the coarse raster) of each cell, from its id (row 0 = south)
function block_of(i)
    v = parse(UInt64, MESH.cells[i].id; base = 16)
    rw = Int(v >> 32); c = Int(v & 0xffffffff)
    return c, (NY ÷ K) - 1 - rw
end
pixels_of(i) = (b = block_of(i); [ZP[K*b[1]+cc+1, K*b[2]+rr+1] for cc in 0:(K-1), rr in 0:(K-1)])

@testset "X1 block means from the raster's origin" begin
    @test NC == (NX ÷ K) * (NY ÷ K)
    zc = ArchGDAL.read(COARSE) do ds
        @test ArchGDAL.getgeotransform(ds) ≈ [X0, K * P, 0.0, Y0, 0.0, -K * P]
        ArchGDAL.read(ds, 1)
    end
    @test size(zc) == (NX ÷ K, NY ÷ K)
    @test all(isapprox(zc[c+1, r+1], sum(ZP[K*c+1:K*c+K, K*r+1:K*r+K]) / K^2; atol = 1e-12)
              for c in 0:(NX ÷ K - 1), r in 0:(NY ÷ K - 1))
    # the mesh's elevation is the block mean
    @test all(isapprox(MESH.static_vars["elevation"][i], sum(pixels_of(i)) / K^2; atol = 1e-9) for i in 1:NC)
end

@testset "X2 cell tables hold exactly what the pixels store" begin
    m = deepcopy(MESH)
    _quiet(() -> build_sgs_tables!(m, FINE; sampling = :pixels, n_bins = 10))
    @test size(m.array_vars["sgs_elev_bins"], 1) == K^2 + 2           # raised to hold every pixel
    @test all(==(1.0), m.static_vars["sgs_sampling"])
    worst = 0.0; worst_eta = 0.0
    rng = MersenneTwister(3)
    for i in 1:NC
        t = sgs_table(m, i); zp = vec(pixels_of(i)); dA = t.cell_area / K^2
        @test t.z_min == minimum(zp)
        @test t.area_curve[end] ≈ t.cell_area rtol = 1e-12
        vexact(η) = sum(dA * max(η - z, 0.0) for z in zp)
        for η in vcat(zp, t.z_min .+ (t.z_max - t.z_min + 1.0) .* rand(rng, 20))
            V = vexact(η)
            V > 0 || continue
            worst = max(worst, abs(A5Grid.wse_from_volume(t, V) - η))
            # the table's own volume at η (linear between bins)
            k = searchsortedlast(t.elev_bins, η)
            Vt = k >= length(t.elev_bins) ? t.vol_curve[end] + (η - t.elev_bins[end]) * t.cell_area :
                 t.vol_curve[k] + (η - t.elev_bins[k]) / (t.elev_bins[k+1] - t.elev_bins[k]) *
                                  (t.vol_curve[k+1] - t.vol_curve[k])
            worst_eta = max(worst_eta, abs(Vt - V) / max(V, 1.0))
        end
    end
    @test worst < 1e-9           # m
    @test worst_eta < 1e-12      # relative volume
    # bilinear sampling (the default) does not: it describes a smoother surface
    mb = deepcopy(MESH)
    _quiet(() -> build_sgs_tables!(mb, FINE; n_samples = 400))
    @test all(==(0.0), mb.static_vars["sgs_sampling"])
    t = sgs_table(mb, 8); zp = vec(pixels_of(8))
    @test t.z_min > minimum(zp) + 0.01
end

@testset "X3 face sills: the higher pixel of each pair across the face" begin
    m = deepcopy(MESH)
    _quiet(() -> build_sgs_tables!(m, FINE; sampling = :pixels, connectivity = false))
    sills = m.array_vars["sgs_edge_sills"]
    idof = Dict(c.id => i for (i, c) in enumerate(m.cells))
    nchecked = 0; nbound = 0
    for i in 1:NC
        ci, ri = block_of(i)
        nbrs = m.adjacency[m.cells[i].id]
        for (slot, nb) in enumerate(nbrs)
            j = idof[nb]; cj, rj = block_of(j)
            pairs = if cj == ci            # north/south neighbour: a row boundary
                rin = rj > ri ? K * ri + K - 1 : K * ri; rout = rj > ri ? rin + 1 : rin - 1
                [(K * ci + q, rin, K * ci + q, rout) for q in 0:(K-1)]
            else
                cin = cj > ci ? K * ci + K - 1 : K * ci; cout = cj > ci ? cin + 1 : cin - 1
                [(cin, K * ri + q, cout, K * ri + q) for q in 0:(K-1)]
            end
            want = minimum(max(ZP[a+1, b+1], ZP[c+1, d+1]) for (a, b, c, d) in pairs)
            @test sills[slot, i] ≈ want atol = 1e-12
            # the neighbour holds the same face
            slot_j = findfirst(==(m.cells[i].id), m.adjacency[nb])
            @test sills[slot_j, j] == sills[slot, i]
            nchecked += 1
        end
        # boundary sides: the lowest inside pixel (the west edge of the first column)
        if ci == 0
            free = A5Grid.boundary_side_slots(size(sills, 1), trues(length(nbrs)))
            ins = [ZP[1, K*ri+q+1] for q in 0:(K-1)]
            @test any(s -> isapprox(sills[s, i], minimum(ins); atol = 1e-12), free)
            nbound += 1
        end
    end
    @test nchecked == 2 * ((NX ÷ K - 1) * (NY ÷ K) + (NX ÷ K) * (NY ÷ K - 1))
    @test nbound == NY ÷ K
end

@testset "X5 a flat block (every pixel at one level)" begin
    zf = copy(ZP); zf[1:K, 1:K] .= 12.5                 # the north-west block
    flat = write_tif(joinpath(DIR, "flat.tif"), zf, X0, Y0, P)
    m = deepcopy(MESH)
    _quiet(() -> build_sgs_tables!(m, flat; sampling = :pixels))
    i = findfirst(j -> block_of(j) == (0, 0), 1:NC)
    t = sgs_table(m, i)
    @test t.z_min == 12.5 && all(==(12.5), t.elev_bins)
    @test A5Grid.wse_from_volume(t, 2.0 * t.cell_area) ≈ 14.5
end

@testset "X4 a mesh not aligned with the DEM is refused" begin
    shifted = write_tif(joinpath(DIR, "shifted.tif"), ZP, X0 + 2.5, Y0, P)
    m = deepcopy(MESH)
    @test_throws ErrorException _quiet(() -> build_sgs_tables!(m, shifted; sampling = :pixels))
    @test_throws ArgumentError build_sgs_tables!(deepcopy(MESH), FINE; sampling = :nearest)
end
