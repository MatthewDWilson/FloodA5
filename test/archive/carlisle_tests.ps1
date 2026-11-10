
# Mesh generation
# Standard res 16 A5 mesh (geometry + neighbours + elevation + cell area)
julia --threads auto --project=. FloodModel.jl `
                --meshgen test/carlisle/Carlisle_domain.geojson `
                --meshres 16 `
                --meshout test/carlisle/carlisle_mesh16_standard.parquet `
                --dem test/carlisle/Carlisle_LiDAR_5m_mean.tif `
                --mesh-only
# SGS res 16 A5 mesh (geometry + neighbours + elevation + cell area)
julia --threads auto --project=. FloodModel.jl `
                --meshgen test/carlisle/Carlisle_domain.geojson `
                --meshres 16 `
                --meshout test/carlisle/carlisle_mesh16_sgs.parquet `
                --dem test/carlisle/Carlisle_LiDAR_5m_mean.tif `
                --mesh-only `
				--flow-model sgs
# Standard res 18 A5 mesh (geometry + neighbours + elevation + cell area)
julia --threads auto --project=. FloodModel.jl `
                --meshgen test/carlisle/Carlisle_domain.geojson `
                --meshres 18 `
                --meshout test/carlisle/carlisle_mesh18_standard.parquet `
                --dem test/carlisle/Carlisle_LiDAR_5m_mean.tif `
                --mesh-only
# SGS res 18 A5 mesh (geometry + neighbours + elevation + cell area)
julia --threads auto --project=. FloodModel.jl `
                --meshgen test/carlisle/Carlisle_domain.geojson `
                --meshres 18 `
                --meshout test/carlisle/carlisle_mesh18_sgs.parquet `
                --dem test/carlisle/Carlisle_LiDAR_5m_mean.tif `
                --mesh-only `
				--flow-model sgs
				
				
# Run test simulations: mesh 16, uniform rainfall, sgs and standard, corrected and uncorrected					

# Run standard simulation res 16 - uniform rainfall (mesh must exist) UNCORRECTED
julia --threads auto --project=. FloodModel.jl `
                --meshload test/carlisle/carlisle_mesh16_standard.parquet `
                --rainfall 50 `
                --output test/carlisle/carlisle_mesh16_standard_unifrom_uncorrected.h5 `
                --output-interval 3600 `
                --sim-duration 36000 `
                --flow-model standard `
				--gradient-correction off *> test/carlisle/mesh16_uniform.txt

# Run sgs simulation res 16 - uniform rainfall (mesh must exist) UNCORRECTED
julia --threads auto --project=. FloodModel.jl `
                --meshload test/carlisle/carlisle_mesh16_sgs.parquet `
                --rainfall 50 `
                --output test/carlisle/carlisle_mesh16_sgs_unifrom_uncorrected.h5 `
                --output-interval 3600 `
                --sim-duration 36000 `
                --flow-model sgs `
				--gradient-correction off *>> test/carlisle/mesh16_uniform.txt

# Run standard simulation res 16 - uniform rainfall (mesh must exist) flow direction corrected
julia --threads auto --project=. FloodModel.jl `
                --meshload test/carlisle/carlisle_mesh16_standard.parquet `
                --rainfall 50 `
                --output test/carlisle/carlisle_mesh16_standard_unifrom_corrected.h5 `
                --output-interval 3600 `
                --sim-duration 36000 `
                --flow-model standard `
				--gradient-correction on *>> test/carlisle/mesh16_uniform.txt

# Run sgs simulation res 16 - uniform rainfall (mesh must exist) flow direction corrected
julia --threads auto --project=. FloodModel.jl `
                --meshload test/carlisle/carlisle_mesh16_sgs.parquet `
                --rainfall 50 `
                --output test/carlisle/carlisle_mesh16_sgs_unifrom_corrected.h5 `
                --output-interval 3600 `
                --sim-duration 36000 `
                --flow-model sgs `
				--gradient-correction on *>> test/carlisle/mesh16_uniform.txt


# Run test simulations: mesh 18, point source, sgs and standard, corrected and uncorrected					

# Run standard simulation res 18 - point source (mesh must exist) UNCORRECTED
julia --threads auto --project=. FloodModel.jl `
                --meshload test/carlisle/carlisle_mesh18_standard.parquet `
				--rainpoint 54.908,-2.896,50.0 `
                --output test/carlisle/carlisle_mesh18_standard_pointsource_uncorrected.h5 `
                --output-interval 3600 `
                --sim-duration 36000 `
                --flow-model standard `
				--gradient-correction off *> test/carlisle/mesh18_pointsource.txt

# Run sgs simulation res 18 - point source (mesh must exist) UNCORRECTED
julia --threads auto --project=. FloodModel.jl `
                --meshload test/carlisle/carlisle_mesh18_sgs.parquet `
				--rainpoint 54.908,-2.896,50.0 `
                --output test/carlisle/carlisle_mesh18_sgs_pointsource_uncorrected.h5 `
                --output-interval 3600 `
                --sim-duration 36000 `
                --flow-model sgs `
				--gradient-correction off *>> test/carlisle/mesh18_pointsource.txt

# Run standard simulation res 18 - point source (mesh must exist) flow direction corrected
julia --threads auto --project=. FloodModel.jl `
                --meshload test/carlisle/carlisle_mesh18_standard.parquet `
				--rainpoint 54.908,-2.896,50.0 `
                --output test/carlisle/carlisle_mesh18_standard_pointsource_corrected.h5 `
                --output-interval 3600 `
                --sim-duration 36000 `
                --flow-model standard `
				--gradient-correction on *>> test/carlisle/mesh18_pointsource.txt

# Run sgs simulation res 18 - point source (mesh must exist) flow direction corrected
julia --threads auto --project=. FloodModel.jl `
                --meshload test/carlisle/carlisle_mesh18_sgs.parquet `
				--rainpoint 54.908,-2.896,50.0 `
                --output test/carlisle/carlisle_mesh18_sgs_pointsource_corrected.h5 `
                --output-interval 3600 `
                --sim-duration 36000 `
                --flow-model sgs `
				--gradient-correction on *>> test/carlisle/mesh18_pointsource.txt

	

 
 