[ Info: A5Grid: CUDA GPU detected ÔÇö PIP sampling will run on GPU
Ôöî Info: Starting FloodA5 model...
Ôöö   Dates.now() = 2026-08-22T19:25:19.597
Ôöî Info: FloodA5 model finished.
Ôöö   Dates.now() = 2026-08-22T19:25:20.912
============================================================================
test_mirror_symmetry.jl
Diagnostic: does the corrected kernel preserve north/south symmetry
on a hand-built, exactly-mirror-symmetric synthetic mesh?
============================================================================

Input WSE (function of longitude only, by construction):
  C: wse=5.0
  E: wse=4.6
  NE: wse=4.85
  SE: wse=4.85
  NW: wse=5.35
  SW: wse=5.35

Edge geometry (by construction, NE/SE and NW/SW are exact mirrors):
  E: dx=400.0 dy=0.0 L=400.0 cos_theta=1.0 skew=(0.0,0.0) nf=(1.0,0.0) sill=3.5999999999999996
  NE: dx=150.0 dy=380.0 L=408.534 cos_theta=0.85 skew=(0.15,0.3) nf=(0.4621,1.0906) sill=3.8499999999999996
  SE: dx=150.0 dy=-380.0 L=408.534 cos_theta=0.85 skew=(0.15,-0.3) nf=(0.4621,-1.0906) sill=3.8499999999999996
  NW: dx=-350.0 dy=300.0 L=460.977 cos_theta=0.8 skew=(-0.1,0.25) nf=(-0.7074,0.7706) sill=4.0
  SW: dx=-350.0 dy=-300.0 L=460.977 cos_theta=0.8 skew=(-0.1,-0.25) nf=(-0.7074,-0.7706) sill=4.0

--- Step A: _build_wlsq_weights! / _compute_wse_gradients! ---
  grad_wse[:,C] = (-0.0009999999999999998, 0.0)
  (gy should be ~0: the input WSE field has zero y-dependence and
   the stencil geometry is exactly mirror-symmetric about y=0)
Test Summary:                                                  | Pass  Total  Time
MS1 -- WLSQ gradient has zero y-component on symmetric stencil |    1      1  0.6s

--- Step B: dWSE_n per edge (Phase A formula, replicated) ---
  E: dWSE_n = 0.40000000000000036
  NE: dWSE_n = 0.15814004732372355
  SE: dWSE_n = 0.15814004732372355
  NW: dWSE_n = -0.303048861143232
  SW: dWSE_n = -0.303048861143232
Test Summary:                                       | Pass  Total  Time
MS2 -- dWSE_n identical for mirror-image edge pairs |    2      2  0.0s

--- Step C: _bates_flux_limited_corrected per edge ---

  dt = 10.0 s:
    E: Q=-4.8069  q_stored=-0.13734
    NE: Q=-1.310092  q_stored=-0.04367
    SE: Q=-1.310092  q_stored=-0.04367
    NW: Q=2.437777  q_stored=0.087063
    SW: Q=2.437777  q_stored=0.087063
Test Summary:                                                    | Pass  Total  Time
MS3 (dt=10.0) -- flux magnitude identical for mirror-image pairs |    2      2  0.0s
    Phase F reconstruction at C: qvec_u=0.124832  qvec_v=0.0
Test Summary:                                     | Pass  Total  Time
MS4 (dt=10.0) -- Phase F qvec_v ~0 at centre cell |    1      1  0.0s

  dt = 0.5 s:
    E: Q=-0.240345  q_stored=-0.006867
    NE: Q=-0.065505  q_stored=-0.002183
    SE: Q=-0.065505  q_stored=-0.002183
    NW: Q=0.121889  q_stored=0.004353
    SW: Q=0.121889  q_stored=0.004353
Test Summary:                                                   | Pass  Total  Time
MS3 (dt=0.5) -- flux magnitude identical for mirror-image pairs |    2      2  0.0s
    Phase F reconstruction at C: qvec_u=0.006242  qvec_v=0.0
Test Summary:                                    | Pass  Total  Time
MS4 (dt=0.5) -- Phase F qvec_v ~0 at centre cell |    1      1  0.0s

============================================================================
If MS1/MS2/MS3/MS4 all PASS at both dt values: the new local-anchor
formula is symmetric in isolation on this synthetic mesh, as
expected -- this is a necessary but not sufficient condition. The
real test is the res-18/res-16 planar-symmetry dt sweep on the
actual simulation: does removing the purely-gradient-driven design
(the missing local anchor identified 2026-08-18) fix the real-mesh
north/south instability, or at least remove its V-shaped, worse-at-
small-dt tail?

If any assertion FAILS: something is still wrong with the new
formula's implementation even in this trivial symmetric case -- the
printed values above pinpoint exactly which term diverges.
============================================================================
