[ Info: A5Grid: CUDA GPU detected ÔÇö PIP sampling will run on GPU
Ôöî Info: Starting FloodA5 model...
Ôöö   Dates.now() = 2026-08-22T19:25:41.672
Ôöî Info: FloodA5 model finished.
Ôöö   Dates.now() = 2026-08-22T19:25:43.028
============================================================
test_cell_momentum.jl
Cell-vector discharge: build, reconstruct, integrate
============================================================

ÔöÇÔöÇ CM1: _build_cell_edge_index! ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
Test Summary:         | Pass  Total  Time
CM1 ÔÇö cell_edge_index |    7      7  1.6s

ÔöÇÔöÇ CM2: _build_mom_weights! ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
Test Summary:                       | Pass  Total  Time
CM2 ÔÇö mom_weights WLSQ construction |    4      4  0.1s

ÔöÇÔöÇ CM2b: _build_mom_weights! ÔÇö regular pentagon ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
   5 face normals at 0┬░, 72┬░, 144┬░, 216┬░, 288┬░;
   inject qvec=(1,0) eastward; verify recovered to < 1% error
   qvec input: (1.000, 0.000) eastward
   qvec recovered: (1.0000, 0.0000)
   qvec input: (0.000, 1.000) northward
   qvec recovered: (0.0000, 1.0000)
Test Summary:                                         | Pass  Total  Time
CM2b ÔÇö WLSQ recovers known vector on regular pentagon |    4      4  0.4s

ÔöÇÔöÇ CM3: Phase F reconstruction from known fluxes ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
   Regular pentagon, flux = 1.0 on eastward face, 0 elsewhere
   Expected: qvec approximately eastward
   flux[1(east)] = -1.0 ÔåÆ qvec = (0.4000, 0.0000)
   flux[4(216┬░)] = -1.0 ÔåÆ qvec = (-0.3236, -0.2351)
Test Summary:                                   | Pass  Total  Time
CM3 ÔÇö Phase F reconstruction sign and direction |    3      3  0.2s

ÔöÇÔöÇ CM4: step_standard! integration with momentum_model=:cell ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
   3-cell chain, cell 1 has higher WSE ÔåÆ should lose volume
   qvec should be non-zero after step
[ Info: step_standard! call 1: NaN_vol=0  NaN_elev=0 NaN_area=0  zero_area=0 vol_sum=40000.0  max_water_depth=3.0
[ Info:   cell[1]: vol=30000.0  depth=3.0  wse=5.0
[ Info:   cell[2]: vol=10000.0  depth=1.0  wse=2.0
[ Info:   cell[3]: vol=0.0  depth=0.0  wse=0.0
   vol[1] 30000.0ÔåÆ27000.0  vol[2] 10000.0ÔåÆ12000.0  vol[3] 0.0ÔåÆ1000.0
   qvec_u = [0.0000, 0.0000, 0.0000]
   qvec_v = [0.0000, 0.0000, 0.0000]
Test Summary:                                | Pass  Total  Time
CM4 ÔÇö step_standard! cell-vector integration |    7      7  1.2s

ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
Results: all @testset assertions passed Ô£ô
