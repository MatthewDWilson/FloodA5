# test/wp1_frozen_isotropy_mesh.jl
# --------------------------------
# WP1 diagnostic: frozen-coefficient friction-limit isotropy on a real mesh
# (A5 or Cartesian GeoParquet), for every combination of face gradient
# {legacy two-point, diamond} × friction {component, vector} × momentum
# {edge, cell}. Reproduces the review's §3.2 table in FloodA5 itself.
# See test/frozen_isotropy_lib.jl for the method.
#
#   FLOODA5_NO_VIS=1 julia --project=. -t 8 test/wp1_frozen_isotropy_mesh.jl \
#       test/planar_embankment/planar_mesh16_std.parquet [more.parquet ...]

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
include(joinpath(@__DIR__, "frozen_isotropy_lib.jl"))

using Printf, Logging

# "edge" uses Q-centring at the default θ = 0.9; "edge1" is θ = 1 (none).
# On pentagons Q-centring mixes non-parallel "collinear" faces, which by
# itself deflects the steady flow (WP1 finding; motivates WP2's θ-blend).
const CONFIGS = [
    ("legacy",  "component", "edge",  (; friction_model = :component)),
    ("legacy",  "vector",    "edge",  (; friction_model = :vector)),
    ("legacy",  "vector",    "edge1", (; friction_model = :vector, q_centre_theta = 1.0)),
    ("diamond", "component", "edge",  (; friction_model = :component, gradient_correction = true, face_flux_method = :diamond)),
    ("diamond", "vector",    "edge",  (; friction_model = :vector,    gradient_correction = true, face_flux_method = :diamond)),
    ("diamond", "component", "edge1", (; friction_model = :component, gradient_correction = true, face_flux_method = :diamond, q_centre_theta = 1.0)),
    ("diamond", "vector",    "edge1", (; friction_model = :vector,    gradient_correction = true, face_flux_method = :diamond, q_centre_theta = 1.0)),
    ("diamond", "component", "cell", (; friction_model = :component, gradient_correction = true, face_flux_method = :diamond, momentum_model = :cell)),
    ("diamond", "vector",    "cell", (; friction_model = :vector,    gradient_correction = true, face_flux_method = :diamond, momentum_model = :cell)),
]

for path in ARGS
    mesh = load_mesh_geoparquet(path)
    println("\n", basename(path), ": ", length(mesh.cells), " cells")
    @printf("  %-8s %-10s %-5s %10s %10s %10s %8s\n", "gradient", "friction", "mom",
            "max|mean|", "rms mean", "cell rms", "cond")
    for (g, f, m, kw) in CONFIGS
        r = with_logger(SimpleLogger(stderr, Logging.Warn)) do
            # dt = 5 s keeps the volume limiter (h·w/(5 dt)) clear of the
            # steady discharge on ~13 m faces; 1200 steps ≈ 30 relaxation times.
            frozen_isotropy(mesh; angles_deg = 0:5:355, dt = 5.0, iters = 1200, kw...)
        end
        @printf("  %-8s %-10s %-5s %9.2f° %9.2f° %9.2f° %8.3f   (interior %d)\n", g, f, m,
                r.max_abs_mean_defl, sqrt(mean(r.mean_defl .^ 2)), mean(r.rms_defl),
                r.cond_ratio, r.n_interior)
    end
end
