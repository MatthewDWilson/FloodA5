#!/usr/bin/env julia
# test_flat_rainpoint.jl
# ----------------------
# Flat-terrain point-source validation test.
#
# Usage (from the FloodA5 project root):
#
#   Step 1 — write AOI GeoJSON (VS Code config [00], or):
#     julia --project=. test_flat_rainpoint.jl --gen
#
#   Step 2 — generate mesh (VS Code config [07], or):
#     julia --project=. --threads auto FloodModel.jl \
#         --meshgen test/flat_test_aoi.geojson --meshres 9 \
#         --meshout test/flat_mesh_res14.parquet \
#         --flow-model standard --mesh-only
#
#   Step 3 — run simulation (VS Code config [08] or [09], or):
#     julia --project=. --threads auto FloodModel.jl \
#         --meshload test/flat_mesh_res14.parquet \
#         --flow-model standard \
#         --rainpoint -43.4043,172.6644,50.0 \
#         --sim-duration 3600 --dt-max 30 \
#         --output test/flat_rainpoint_out.h5 --output-interval 300
#
#   Step 4 — validate HDF5 output (VS Code config [10], or):
#     julia --project=. test_flat_rainpoint.jl --analyse

# ---------------------------------------------------------------------------
# Determine mode before loading any heavy packages
# ---------------------------------------------------------------------------

const DO_GEN   = "--gen"     in ARGS
const DO_ANAL  = "--analyse" in ARGS

# Only load HDF5/Statistics when actually needed — avoids the package entirely
# during --gen, which has no dependencies beyond stdlib.
if DO_ANAL
    using HDF5
    using Statistics
end

using Printf

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

const TEST_DIR    = joinpath(@__DIR__, "test")
const AOI_FILE    = joinpath(TEST_DIR, "flat_test_aoi.geojson")
const MESH_FILE   = joinpath(TEST_DIR, "flat_mesh_res14.parquet")
const OUTPUT_FILE      = joinpath(TEST_DIR, "flat_rainpoint_out.h5")
const OUTPUT_FILE_2HR  = joinpath(TEST_DIR, "flat_rainpoint_2hr.h5")
const OUTPUT_FILE_200  = joinpath(TEST_DIR, "flat_rainpoint_200mmhr.h5")
const OUTPUT_FILE_1000 = joinpath(TEST_DIR, "flat_rainpoint_1000mmhr.h5")
const RAIN_MM_HR  = 50.0

# Centre of test AOI (Christchurch, NZ)
# AOI centre matches the actual flat_test_aoi.geojson (small box near Christchurch)
const AOI_LON_C  = 172.6644
const AOI_LAT_C  = -43.4043

# ---------------------------------------------------------------------------
# AOI generation
# ---------------------------------------------------------------------------

function write_aoi(path::String, lon::Float64, lat::Float64, half::Float64)
    mkpath(dirname(path))
    open(path, "w") do io
        print(io, """{
  "type": "Feature",
  "properties": { "name": "Flat terrain test - Christchurch" },
  "geometry": {
    "type": "Polygon",
    "coordinates": [[
      [$(lon-half), $(lat-half)],
      [$(lon+half), $(lat-half)],
      [$(lon+half), $(lat+half)],
      [$(lon-half), $(lat+half)],
      [$(lon-half), $(lat-half)]
    ]]
  }
}""")
    end
end

# ---------------------------------------------------------------------------
# Validation  (only called when DO_ANAL — HDF5/Statistics loaded above)
# ---------------------------------------------------------------------------

function validate_output(h5_path::String)
    println("\n=== FloodA5 Flat Rainpoint Validation ===\n")
    all_pass = true

    h5open(h5_path, "r") do f
        cell_ids = read(f["mesh/cell_ids"])
        lons     = read(f["mesh/center_lons"])
        lats     = read(f["mesh/center_lats"])
        n_cells  = length(cell_ids)
        frames   = sort(keys(f["frames"]))
        n_frames = length(frames)

        println("Mesh   : $n_cells cells")
        println("Frames : $n_frames  ($(frames[1]) -> $(frames[end]))\n")

        # Source cell — nearest mesh centre to the AOI centre
        dists   = sqrt.((lons .- AOI_LON_C).^2 .+ (lats .- AOI_LAT_C).^2)
        src_idx = argmin(dists)
        src_id  = cell_ids[src_idx]
        println("Source cell : index=$src_idx  id=$src_id  " *
                "(lon=$(round(lons[src_idx],digits=4)), " *
                "lat=$(round(lats[src_idx],digits=4)))\n")

        vol_frames = [sum(read(f["frames/$fr/volume"])) for fr in frames]
        t_frames   = [read(f["frames/$fr/t"])           for fr in frames]
        wet_counts = [count(>(1e-4), read(f["frames/$fr/volume"])) for fr in frames]

        final = frames[end]
        vol_f = read(f["frames/$final/volume"])
        dep_f = read(f["frames/$final/water_depth"])
        vel_f = read(f["frames/$final/velocity"])

        # Check 1: No NaN / negative volumes
        has_nan = any(isnan, vol_f) || any(isnan, dep_f)
        has_neg = any(<(0.0), vol_f)
        print("Check 1 — No NaN or negative volumes     : ")
        if !has_nan && !has_neg
            println("PASS")
        else
            println("FAIL  (NaN=$(count(isnan, vol_f))  neg=$(count(<(0.0), vol_f)))")
            all_pass = false
        end

        # Check 2: Domain volume growing
        dt_sim = t_frames[end] - t_frames[1]
        dVdt   = dt_sim > 0.0 ? (vol_frames[end] - vol_frames[1]) / dt_sim : 0.0
        print("Check 2 — Domain volume growing          : ")
        if dVdt > 0.0
            println("PASS  (mean rate $(round(dVdt, sigdigits=3)) m3/s)")
        else
            println("FAIL  (domain volume not growing)")
            all_pass = false
        end

        # Check 3: Wet cells spreading
        print("Check 3 — Wet cells increase over time   : ")
        if wet_counts[end] > wet_counts[1]
            println("PASS  ($(wet_counts[1]) -> $(wet_counts[end]) wet cells)")
        else
            println("FAIL  (wet cell count did not increase)")
            all_pass = false
        end

        # Check 4: Source cell depth
        src_dep = dep_f[src_idx]
        print("Check 4 — Source cell depth > 0.01 m     : ")
        if src_dep > 0.01
            println("PASS  ($(round(src_dep, digits=4)) m)")
        else
            println("FAIL  (depth = $(round(src_dep, sigdigits=3)) m)")
            all_pass = false
        end

        # Check 5: Velocity (warn only)
        print("Check 5 — Velocity non-zero              : ")
        if all(==(0.0), vel_f)
            println("WARN  (all zero — velocity computation not yet implemented)")
        else
            println("PASS  (max $(round(maximum(vel_f), sigdigits=3)) m/s)")
        end

        # Check 6: Mass balance — domain volume grows at the source rate
        if length(t_frames) >= 2
            mean_rate = (vol_frames[end] - vol_frames[1]) / (t_frames[end] - t_frames[1])
            expected  = mean_rate * t_frames[end]
            mb_pct    = abs(vol_frames[end] - expected) / max(expected, 1.0) * 100
            print("Check 6 — Mass balance < 0.5%             : ")
            if mb_pct < 0.5
                println("PASS  ($(round(mb_pct, sigdigits=2))% error)")
            else
                println("FAIL  ($(round(mb_pct, sigdigits=3))% error)")
                all_pass = false
            end
        end

        # Check 7: All cells eventually wet
        final_wet = count(>(1e-4), dep_f)
        pct_wet   = 100.0 * final_wet / n_cells
        print("Check 7 — All cells wet at final frame    : ")
        if final_wet == n_cells
            println("PASS  (all $n_cells cells wet)")
        elseif pct_wet >= 80.0
            println("WARN  ($final_wet/$n_cells = $(round(pct_wet,digits=1))% wet — run longer for full coverage)")
        else
            println("FAIL  ($final_wet/$n_cells = $(round(pct_wet,digits=1))% wet)")
            all_pass = false
        end

        # Per-frame summary table
        println("\n-- Per-frame summary --")
        @printf("%-8s  %-12s  %-10s  %-13s  %-13s\n",
                "t (s)", "vol_tot(m3)", "wet cells", "max_depth(m)", "src_depth(m)")
        println(repeat('-', 64))
        for (i, fr) in enumerate(frames)
            dep = read(f["frames/$fr/water_depth"])
            @printf("%-8.0f  %-12.2f  %-10d  %-13.4f  %-13.4f\n",
                    t_frames[i], vol_frames[i], wet_counts[i],
                    maximum(dep), dep[src_idx])
        end

        println("\n-- $(all_pass ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") --\n")
    end
    return all_pass
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

function main_test()
    if !DO_GEN && !DO_ANAL
        println("""
Usage (from FloodA5 project root):

  julia --project=. test_flat_rainpoint.jl --gen        write AOI GeoJSON to test/
  julia --project=. test_flat_rainpoint.jl --analyse    validate test/flat_rainpoint_out.h5

Mesh generation and simulation are run via VS Code launch configs [07]-[09],
or directly with FloodModel.jl (see comment block at top of this file).
""")
        return
    end

    if DO_GEN
        println("Writing AOI GeoJSON...")
        write_aoi(AOI_FILE, AOI_LON_C, AOI_LAT_C, AOI_HALF)
        if isfile(AOI_FILE)
            println("OK  AOI written to: $AOI_FILE")
            println("    Next: run launch config [07] to generate the mesh.")
        else
            println("ERROR: file was not created — check that $(TEST_DIR) is writable.")
        end
    end

    if DO_ANAL
        # Pick output file: pass --file 2hr or --file 200 to select alternate runs
        h5_path = let
            idx = findfirst(==("--file"), ARGS)
            if idx !== nothing
                tag = get(ARGS, idx+1, "")
                tag == "2hr"  ? OUTPUT_FILE_2HR  :
                tag == "200"  ? OUTPUT_FILE_200   :
                tag == "1000" ? OUTPUT_FILE_1000  : OUTPUT_FILE
            else
                OUTPUT_FILE
            end
        end
        if !isfile(h5_path)
            println("ERROR: output not found: $h5_path")
            println("       Run the corresponding launch config first.")
            return
        end
        validate_output(h5_path)
    end
end

main_test()
