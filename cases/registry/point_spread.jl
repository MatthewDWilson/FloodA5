# cases/registry/point_spread.jl — flat-bed point source in a closed square
# (docs/TEST_CASES.md §1, §5b). The spreading pattern should be circular; its
# departures measure the mesh's directional bias (A4: four-fold, Cartesian
# diagonal vs cardinal; A6: six-fold; R: cardinal/diagonal front radius).

const _SQUARE_AOI = "{repo}/test/square/square_domain.geojson"   # 1.5 km square at (0, 0)

"Point-spread metrics at the last frame of every run."
function _point_spread_metrics(ctx)
    out = Dict{String,Dict{String,Float64}}()
    for (name, r) in ctx.runs
        k = FA.nframes(r)
        s = FA.spread_metrics(r, k; lon0 = 0.0, lat0 = 0.0)
        p = FA.point_spread_shape(r, k; lon0 = 0.0, lat0 = 0.0, saturation_threshold = 1.0)
        out[name] = Dict("t" => s.t, "n_wet" => s.n_wet, "A4" => s.A4, "A6" => s.A6,
                         "ratio_1cm" => s.ratio[2], "circularity" => p.circularity,
                         "cv_radius" => p.cv_radius, "storage_m3" => FA.storage(r, k))
    end
    return out
end

_pt_all(f, m) = all(f(d) for d in values(m))
_pt_worst(m, key) = join(("$j $(round(d[key]; sigdigits = 3))" for (j, d) in sort(collect(m); by = first)), ", ")

register!(Case(
    name = "point_spread_smoke", tier = :smoke,
    title = "Flat point source, Cartesian 50 m, 20 min (CI)",
    doc = "docs/TEST_CASES.md §5b",
    description = """
    The smallest end-to-end case: a 900-cell Cartesian grid on the square AOI,
    2 m³/s injected at the centre for 20 minutes with closed boundaries. Runs
    in CI through the runner (prepare, run, analyse, check against expected).""",
    prep = [Prep("mesh", :floodmodel,
                 ["--meshgen-cartesian", _SQUARE_AOI, "--dx", "50",
                  "--meshout", "{case}/mesh_dx50.parquet", "--mesh-only"];
                 outputs = ["{case}/mesh_dx50.parquet"])],
    jobs = [Job("cart_dx50", ["--meshload", "{case}/mesh_dx50.parquet", "--flow-model", "standard",
                               "--injection-point", "0.0,0.0,2.0", "--closed-boundaries",
                               "--sim-duration", "1200", "--dt-max", "10", "--output-interval", "300"];
                threads = 1, mesh = "cartesian", dx = 50)],
    analyse = _point_spread_metrics,
    checks = [Check("mass conserved (|error| < 1e-9 of the volume added)",
                    m -> _pt_all(d -> d["mass_balance_rel"] < 1e-9, m); detail = m -> _pt_worst(m, "mass_balance_rel")),
              Check("stored volume = 2 m³/s × 1200 s (±0.1%)",
                    m -> _pt_all(d -> abs(d["storage_m3"] / 2400 - 1) < 1e-3, m); detail = m -> _pt_worst(m, "storage_m3"))],
    # CI runs on Linux; the expected values are recorded on Windows
    rtol = Dict("*" => 1e-4), atol = Dict("*" => 1e-9),
))

let jobs = Job[]
    for (mesh, res) in (("a5", 16), ("a5", 18))
        push!(jobs, Job("a5_res$res", ["--meshload", "{case}/a5_res$res.parquet", "--flow-model", "standard",
                                       "--injection-point", "0.0,0.0,10.0", "--closed-boundaries",
                                       "--sim-duration", "3600", "--dt-max", "2", "--output-interval", "600"];
                        mesh = mesh, res = res))
    end
    for dx in ("25", "12.5")
        tag = replace(dx, "." => "p")
        push!(jobs, Job("cart_dx$tag", ["--meshload", "{case}/cart_dx$tag.parquet", "--flow-model", "standard",
                                        "--injection-point", "0.0,0.0,10.0", "--closed-boundaries",
                                        "--sim-duration", "3600", "--dt-max", "2", "--output-interval", "600"];
                        mesh = "cartesian", dx = parse(Float64, dx)))
    end
    register!(Case(
        name = "point_spread", tier = :standard,
        title = "Flat point source on A5 and Cartesian meshes (directional bias)",
        doc = "docs/TEST_CASES.md §1, §5b",
        description = """
        10 m³/s injected at the centre of the closed, flat 1.5 km square for one
        hour, standard flow with the default scheme, on A5 (res 16, 18) and
        Cartesian (25, 12.5 m) meshes. Metrics at t = 3600 s: A4, A6 (depth-
        weighted angular harmonics about the source; 0 for a circle), the
        cardinal/diagonal front-radius ratio at 1 cm, hull circularity and the
        coefficient of variation of the front radius over 8 bearings.""",
        prep = vcat([Prep("a5_res$res", :floodmodel,
                          ["--meshgen", _SQUARE_AOI, "--meshres", "$res",
                           "--meshout", "{case}/a5_res$res.parquet", "--flow-model", "standard", "--mesh-only"];
                          outputs = ["{case}/a5_res$res.parquet"]) for res in (16, 18)],
                    [Prep("cart_dx$(replace(dx, "." => "p"))", :floodmodel,
                          ["--meshgen-cartesian", _SQUARE_AOI, "--dx", dx,
                           "--meshout", "{case}/cart_dx$(replace(dx, "." => "p")).parquet", "--mesh-only"];
                          outputs = ["{case}/cart_dx$(replace(dx, "." => "p")).parquet"]) for dx in ("25", "12.5")]),
        jobs = jobs,
        analyse = _point_spread_metrics,
        checks = [Check("mass conserved (|error| < 1e-9 of the volume added)",
                        m -> _pt_all(d -> d["mass_balance_rel"] < 1e-9, m); detail = m -> _pt_worst(m, "mass_balance_rel")),
                  Check("near-circular spread (hull circularity > 0.9)",
                        m -> _pt_all(d -> d["circularity"] > 0.9, m); detail = m -> _pt_worst(m, "circularity")),
                  Check("weak four-fold bias (|A4| < 0.1)",
                        m -> _pt_all(d -> abs(d["A4"]) < 0.1, m); detail = m -> _pt_worst(m, "A4"))],
        rtol = Dict("*" => 1e-6), atol = Dict("*" => 1e-9),
    ))
end
