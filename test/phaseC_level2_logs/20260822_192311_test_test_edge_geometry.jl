[ Info: A5Grid: CUDA GPU detected ÔÇö PIP sampling will run on GPU
Ôöî Info: Starting FloodA5 model...
Ôöö   Dates.now() = 2026-08-22T19:23:27.848
Ôöî Info: FloodA5 model finished.
Ôöö   Dates.now() = 2026-08-22T19:23:27.993

Group 1 ÔÇö _edge_geometry (cos ╬© component) unit geometry
ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
  [32mÔ£ô[0m  1.1  orthogonal pair ÔåÆ cos ╬© Ôëê 1.0
  [32mÔ£ô[0m  1.2  45┬░ skew ÔåÆ cos ╬© Ôëê 0.7071
  [32mÔ£ô[0m  1.3  20┬░ skew ÔåÆ cos ╬© Ôëê 0.9397
  [32mÔ£ô[0m  1.4  non-adjacent cells ÔåÆ 1.0 fallback (no NaN suppression)
  [32mÔ£ô[0m  1.5  large centre separation ÔåÆ finite result, no crash
  [32mÔ£ô[0m  1.6  symmetry: swap iÔåöj ÔåÆ same cos ╬©
  [32mÔ£ô[0m  1.7  resolution invariance: coarse Ôëê fine cos ╬©

Group 2 ÔÇö _build_edge_list and EdgeList population
ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
[ Info: Edge list built: 2 edges for 3 cells (0.67 edges/cell)
[ Info: Edge non-orthogonality (cos ╬©):  min=1.0  mean=1.0  max=1.0
[ Info: Edge non-orthogonality (╬©, degrees):  min=0.0┬░  mean=0.0┬░  p95=0.0┬░  max=0.0┬░  (OpenFOAM: corrected-scheme-safe Ôëñ70┬░, limiting advised above)
[ Info: Edge WLSQ correction vector |V╠é| (= sin ╬©):  mean=0.0  max=0.0  (0 = orthogonal/no correction, bounded in [0,1])
[ Info: Q-centred collinear edges: ci=0/2  cj=0/2 found
  [32mÔ£ô[0m  2.1  edge count correct (3 cells ÔåÆ 2 edges)
  [32mÔ£ô[0m  2.2  all edge widths finite
  [32mÔ£ô[0m  2.2  all edge distances finite
  [32mÔ£ô[0m  2.2  all cos_theta finite and in [0,1]
  [32mÔ£ô[0m  2.2  all sill values finite
  [32mÔ£ô[0m  2.3  canonical ordering: cell_i < cell_j for all edges
  [32mÔ£ô[0m  2.4  centre_dist > 0 for all edges
[ Info: Edge list built: 1 edges for 2 cells (0.5 edges/cell)
[ Info: Edge non-orthogonality (cos ╬©):  min=1.0  mean=1.0  max=1.0
[ Info: Edge non-orthogonality (╬©, degrees):  min=0.0┬░  mean=0.0┬░  p95=0.0┬░  max=0.0┬░  (OpenFOAM: corrected-scheme-safe Ôëñ70┬░, limiting advised above)
[ Info: Edge WLSQ correction vector |V╠é| (= sin ╬©):  mean=0.0  max=0.0  (0 = orthogonal/no correction, bounded in [0,1])
[ Info: Q-centred collinear edges: ci=0/1  cj=0/1 found
  [32mÔ£ô[0m  2.5  mixed-resolution pair ÔåÆ 1 finite edge (MR forward-compat.)

Group 3 ÔÇö _bates_flux integration
ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
  [32mÔ£ô[0m  3.1  cos ╬© = 1.0 reproduces pre-correction formula
  [32mÔ£ô[0m  3.2  cos ╬© < 1 increases |flux| relative to orthogonal
  [32mÔ£ô[0m  3.3  no flow when h_flow Ôëñ 0
  [32mÔ£ô[0m  3.4  pathological skew (cos ╬© = 0.05) ÔåÆ finite flux
  [32mÔ£ô[0m  3.5  sign: Q < 0 when wse_i > wse_j (flow goes iÔåÆj)
  [32mÔ£ô[0m  3.6  sign: Q > 0 when wse_j > wse_i (flow goes jÔåÆi, i gains)

Group 4 ÔÇö EdgeList flux, sign convention, volume update
ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
[ Info: step_standard! call 1: NaN_vol=0  NaN_elev=0 NaN_area=0  zero_area=0 vol_sum=3.75e6  max_water_depth=2.0
[ Info:   cell[1]: vol=3.0e6  depth=2.0  wse=2.0
[ Info:   cell[2]: vol=750000.0  depth=0.5  wse=0.5
  [32mÔ£ô[0m  4.1  single edge: deeper cell loses, shallower cell gains
[ Info: step_standard! call 2: NaN_vol=0  NaN_elev=0 NaN_area=0  zero_area=0 vol_sum=3.75e6  max_water_depth=2.0
[ Info:   cell[1]: vol=3.0e6  depth=2.0  wse=2.0
[ Info:   cell[2]: vol=750000.0  depth=0.5  wse=0.5
  [32mÔ£ô[0m  4.2  NaN cos_theta edge skipped; valid edge still processed
  [32mÔ£ô[0m  4.3  step_standard!: deeper cell loses volume
  [32mÔ£ô[0m  4.4  step_standard!: global volume conserved
  [32mÔ£ô[0m  4.5  step_sgs!: deeper cell loses volume
  [32mÔ£ô[0m  4.6  step_sgs!: global volume conserved

ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
Results: 27 / 27 passed  ÔÇö all tests passed Ô£ô
