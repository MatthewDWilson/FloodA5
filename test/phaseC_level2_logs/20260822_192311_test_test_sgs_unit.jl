[ Info: A5Grid: CUDA GPU detected ÔÇö PIP sampling will run on GPU
Ôöî Info: Starting FloodA5 model...
Ôöö   Dates.now() = 2026-08-22T19:24:37.157
Ôöî Info: FloodA5 model finished.
Ôöö   Dates.now() = 2026-08-22T19:24:37.970

============================================================
SGS unit test ÔÇö 5-cell linear chain
============================================================
  SGSTable lookups: OK
Test Summary:    | Pass  Total  Time
SGSTable lookups |   40     40  1.8s
  Bug 48 refinement: max(z_sill, z_min) dry-cell check ÔÇö OK
    vol after 1 step: [500.0, 0.0, 0.0, 0.0, 0.0]
Test Summary:                                    | Pass  Total  Time
Bug 48 ÔÇö dry cells: wse_eff = max(z_sill, z_min) |    5      5  1.7s
  Mass conservation:
    injected = 1500.0 m┬│
    domain   = 1500.0 m┬│
    error    = 0.0 m┬│  (0.0%)
Test Summary:                                     | Pass  Total  Time
Mass conservation ÔÇö 100 steps, constant injection |    1      1  0.0s
  Volume profile after 300 steps:
    cell 1 (z_min=5.0 m): 1593.65 m┬│
    cell 2 (z_min=6.0 m): 1412.66 m┬│
    cell 3 (z_min=7.0 m): 820.73 m┬│
    cell 4 (z_min=8.0 m): 614.84 m┬│
    cell 5 (z_min=9.0 m): 58.12 m┬│
  Downhill flow: source volume ÔëÑ far cell ÔÇö OK
Test Summary:                              | Pass  Total  Time
Downhill flow ÔÇö source retains most volume |    1      1  0.1s
  No sloshing: max cell volume = 200.0 m┬│ (started with 200.0 m┬│ total) ÔÇö OK
Test Summary:               | Pass  Total  Time
No sloshing ÔÇö free drainage |    2      2  0.0s

============================================================
All SGS unit tests passed.
============================================================

