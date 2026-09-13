julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh18_std.parquet `
    --flow-model standard --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res18_a00.h5 `
    --output-interval 300 `
    --gradient-correction on --gradient-correction-alpha 0.0

julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh18_std.parquet `
    --flow-model standard --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res18_a01.h5 `
    --output-interval 300 `
    --gradient-correction on --gradient-correction-alpha 0.1

julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh18_std.parquet `
    --flow-model standard --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res18_a02.h5 `
    --output-interval 300 `
    --gradient-correction on --gradient-correction-alpha 0.2

julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh18_std.parquet `
    --flow-model standard --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res18_a05.h5 `
    --output-interval 300 `
    --gradient-correction on --gradient-correction-alpha 0.5

julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh18_std.parquet `
    --flow-model standard --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res18_a10.h5 `
    --output-interval 300 `
    --gradient-correction on --gradient-correction-alpha 1.0


    julia --project=. test\test_planar_symmetry.jl `
    --baseline  test/planar_embankment/planar_res18_std_baseline.h5 `
    --corrected test/planar_embankment/planar_res18_a00.h5 `
    --source-lat 51.0001 --source-lon -0.0434 --sweep 10

julia --project=. test\test_planar_symmetry.jl `
    --baseline  test/planar_embankment/planar_res18_std_baseline.h5 `
    --corrected test/planar_embankment/planar_res18_a01.h5 `
    --source-lat 51.0001 --source-lon -0.0434 --sweep 10

julia --project=. test\test_planar_symmetry.jl `
    --baseline  test/planar_embankment/planar_res18_std_baseline.h5 `
    --corrected test/planar_embankment/planar_res18_a02.h5 `
    --source-lat 51.0001 --source-lon -0.0434 --sweep 10

julia --project=. test\test_planar_symmetry.jl `
    --baseline  test/planar_embankment/planar_res18_std_baseline.h5 `
    --corrected test/planar_embankment/planar_res18_a05.h5 `
    --source-lat 51.0001 --source-lon -0.0434 --sweep 10

julia --project=. test\test_planar_symmetry.jl `
    --baseline  test/planar_embankment/planar_res18_std_baseline.h5 `
    --corrected test/planar_embankment/planar_result_std_corrected.h5 `
    --source-lat 51.0001 --source-lon -0.0434 --sweep 10