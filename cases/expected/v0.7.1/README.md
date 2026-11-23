# Expected values as of 0.7.1

The case metrics recorded up to 0.7.1, kept as the record when 0.8.0
re-recorded every case that the DEM registration fix (DEM sampled at pixel
centres, not half a pixel east and south) or the SGS critical-flow limit
(`--sgs-critical`, default on) changes. See `CHANGELOG.md` [0.8.0].

To reproduce one, check out tag `v0.7.1` and run the case there
(`julia --project=. cases/run.jl run CASE`). SGS meshes built by 0.7.x (table
layout 3) are refused by 0.8.0 and must be rebuilt.
