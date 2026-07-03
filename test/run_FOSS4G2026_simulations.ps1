
<#
PowerShell script to run simulations for the FloodA5 paper presented at 
FOSS4G 2026 Hiroshima. These are assumed to be run from the root directory 
of FloodA5 repo, which is available from:
https://github.com/MatthewDWilson/FloodA5/tree/FOSS4G2026
(this is the same version used for the simulations in the paper). 
The libraries, especially the Python bridge and pya5 should be installed 
according to setup.jl. 
#>

#-------------------------------------------------------------------------
# Point spread function: direction uncorrected
julia --threads auto --project=. FloodModel.jl `
	--meshgen test/square/square_domain.geojson `
	--meshres 18 `
	--meshout test/square/square_mesh18_standard.parquet `
    --flow-model standard `
	--rainpoint 0.0,0.0,50.0 `
    --sim-duration 72000 `
	--dt-max 10 `
	--output test/square/square_rainpoint_2hr.h5 `
	--output-interval 300 `
    --gradient-correction off

# Point spread function: direction corrected - updated July 2026
julia --threads auto --project=. FloodModel.jl `
	--meshgen test/square/square_domain.geojson `
	--meshres 18 `
	--meshout test/square/square_mesh18_standard.parquet `
    --flow-model standard `
	--rainpoint 0.0,0.0,50.0 `
    --sim-duration 72000 `
	--dt-max 10 `
	--output test/square/square_rainpoint_2hr_corrected.h5 `
	--output-interval 300 `
    --gradient-correction on

#-------------------------------------------------------------------------
# Planar embankment simulations - gradient correction off (as per FOSS4G2026 paper)
#
# Standard flow model:
julia --project=. --threads auto FloodModel.jl `
              --meshload test/planar_embankment/planar_mesh18_std.parquet `
              --flow-model standard `
              --injection-point 51.0001,-0.0434,0.1 `
              --closed-boundaries `
              --sim-duration 72000 `
              --dt-max 10 `
              --output  test/planar_embankment/planar_result_std.h5 `
              --output-interval 300 `
              --gradient-correction off `
              --vis makie

# SGS flow model:
julia --project=. --threads auto FloodModel.jl `
              --meshload test/planar_embankment/planar_mesh18_sgs.parquet `
              --flow-model sgs `
              --injection-point 51.0001,-0.0434,0.1 `
              --closed-boundaries `
              --sim-duration 72000 `
              --dt-max 10 `
              --output  test/planar_embankment/planar_result_sgs.h5 `
              --output-interval 300 `
              --gradient-correction off `
              --vis makie

# Planar embankment simulations - gradient correction ON (updated July 2026)
#
# Standard flow model:
julia --project=. --threads auto FloodModel.jl `
              --meshload test/planar_embankment/planar_mesh18_std.parquet `
              --flow-model standard `
              --injection-point 51.0001,-0.0434,0.1 `
              --closed-boundaries `
              --sim-duration 72000 `
              --dt-max 10 `
              --output  test/planar_embankment/planar_result_std_corrected.h5 `
              --output-interval 300 `
              --gradient-correction on `
              --vis makie

# SGS flow model:
julia --project=. --threads auto FloodModel.jl `
              --meshload test/planar_embankment/planar_mesh18_sgs.parquet `
              --flow-model sgs `
              --injection-point 51.0001,-0.0434,0.1 `
              --closed-boundaries `
              --sim-duration 72000 `
              --dt-max 10 `
              --output  test/planar_embankment/planar_result_sgs_corrected.h5 `
              --output-interval 300 `
              --gradient-correction on `
              --vis makie


#-------------------------------------------------------------------------
# Carlisle test cases: grid generation (run separetely here but can be
# included as part of a full simulation.

# Mesh resolution 18, standard flow:
julia --threads auto --project=. FloodModel.jl `
	--meshgen test/carlisle/Carlisle_domain.geojson `
	--meshres 18 `
	--meshout test/carlisle/carlisle_mesh18_standard.parquet `
	--dem test/carlisle/Carlisle_LiDAR_5m_mean.tif `
    --flow-model standard `
    --mesh-only
				
# Mesh resolution 18, SGS flow:
julia --threads auto --project=. FloodModel.jl `
	--meshgen test/carlisle/Carlisle_domain.geojson `
	--meshres 18 `
	--meshout test/carlisle/carlisle_mesh18_standard.parquet `
	--dem test/carlisle/Carlisle_LiDAR_5m_mean.tif `
    --flow-model sgs `
    --mesh-only

# Mesh resolution 20, standard flow:
julia --threads auto --project=. FloodModel.jl `
	--meshgen test/carlisle/Carlisle_domain.geojson `
	--meshres 18 `
	--meshout test/carlisle/carlisle_mesh18_standard.parquet `
	--dem test/carlisle/Carlisle_LiDAR_5m_mean.tif `
    --flow-model standard `
    --mesh-only


#-------------------------------------------------------------------------
# Carlisle test cases: model simulations - uncorrected gradient (as per FOSS4G2026 paper)

# Mesh resolution 18, standard flow:
julia --threads auto --project=. FloodModel.jl `
    --meshload test/carlisle/carlisle_mesh18_standard.parquet `
    --flow-model standard `
    --inflow-bci test/carlisle/carlisle.bci `
    --bc-epsg 27700 `
    --manning-n 0.03 `
    --sim-duration 432000 `
    --dt-max 10 `
    --output test/carlisle/carlisle_standard_res18.h5 `
    --output-interval 3600 `
    --gradient-correction off

# Mesh resolution 18, SGS flow:
julia --threads auto --project=. FloodModel.jl `
    --meshload test/carlisle/carlisle_mesh18_sgs.parquet `
    --flow-model sgs `
    --inflow-bci test/carlisle/carlisle.bci `
    --bc-epsg 27700 `
    --manning-n 0.03 `
    --sim-duration 432000 `
    --dt-max 10 `
    --output test/carlisle/carlisle_sgs_res18.h5 `
    --output-interval 3600 `
    --gradient-correction off

# Mesh resolution 20, standard flow:
julia --threads auto --project=. FloodModel.jl `
    --meshload test/carlisle/carlisle_mesh20_standard.parquet `
    --flow-model standard `
    --inflow-bci test/carlisle/carlisle.bci `
    --bc-epsg 27700 `
    --manning-n 0.03 `
    --sim-duration 432000 `
    --dt-max 10 `
    --output test/carlisle/carlisle_standard_res20.h5 `
    --output-interval 3600 `
    --gradient-correction off

# Carlisle test cases: model simulations - corrected gradient (updated July 2026)

# Mesh resolution 18, standard flow:
julia --threads auto --project=. FloodModel.jl `
    --meshload test/carlisle/carlisle_mesh18_standard.parquet `
    --flow-model standard `
    --inflow-bci test/carlisle/carlisle.bci `
    --bc-epsg 27700 `
    --manning-n 0.03 `
    --sim-duration 432000 `
    --dt-max 10 `
    --output test/carlisle/carlisle_standard_res18_corrected.h5 `
    --output-interval 3600 `
    --gradient-correction on

# Mesh resolution 18, SGS flow:
julia --threads auto --project=. FloodModel.jl `
    --meshload test/carlisle/carlisle_mesh18_sgs.parquet `
    --flow-model sgs `
    --inflow-bci test/carlisle/carlisle.bci `
    --bc-epsg 27700 `
    --manning-n 0.03 `
    --sim-duration 432000 `
    --dt-max 10 `
    --output test/carlisle/carlisle_sgs_res18_corrected.h5 `
    --output-interval 3600 `
    --gradient-correction on

# Mesh resolution 20, standard flow:
julia --threads auto --project=. FloodModel.jl `
    --meshload test/carlisle/carlisle_mesh20_standard.parquet `
    --flow-model standard `
    --inflow-bci test/carlisle/carlisle.bci `
    --bc-epsg 27700 `
    --manning-n 0.03 `
    --sim-duration 432000 `
    --dt-max 10 `
    --output test/carlisle/carlisle_standard_res20_corrected.h5 `
    --output-interval 3600 `
    --gradient-correction on
