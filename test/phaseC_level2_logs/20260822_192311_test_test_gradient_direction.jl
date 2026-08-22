[ Info: A5Grid: CUDA GPU detected ÔÇö PIP sampling will run on GPU
==============================================================
test_gradient_direction.jl
d╠é projection formula ÔÇö sign, magnitude, resolution sweep
==============================================================

ÔöÇÔöÇ GD1: V╠é_y for tilted-edge cell pairs, tilt ÔêÆ25┬░ to +25┬░ ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
   Models non-orthogonality seen on real A5 cells (╬© Ôëê 22ÔÇô37┬░)
   Positive tilt (CW from north) ÔåÆ V╠é_y negative; negative tilt ÔåÆ positive
   The sign depends on which way the edge tilts ÔÇö sublattice-dependent on A5

  tilt(┬░)   cos_╬©     V╠é_y       sign_V╠é_y
  ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
  -25.0     0.9063    0.4226    positive
  -20.0     0.9397    0.3420    positive
  -15.0     0.9659    0.2588    positive
  -10.0     0.9848    0.1736    positive
  0.0       1.0000    -0.0000   zero
  10.0      0.9848    -0.1736   negative
  15.0      0.9659    -0.2588   negative
  20.0      0.9397    -0.3420   negative
  25.0      0.9063    -0.4226   negative

  V╠é_y sign counts:  positive=4  negative=4  zero=1
  (On A5 meshes, different cell pairs have different tilt directions,
   giving mixed V╠é_y signs ÔÇö this is why V╠é correction is resolution-dependent)
Test Summary:                                  | Pass  Total  Time
GD1 ÔÇö V╠é_y sign depends on edge tilt direction  |    8      8  0.9s

ÔöÇÔöÇ GD2: d╠é projection formula ÔÇö sign and magnitude, res 12ÔÇô20 ÔöÇÔöÇÔöÇÔöÇ
   WSE field: WSE = ÔêÆ0.001 m/m ├ù x_east  (ci west=uphill, cj east=downhill)
   Expected dWSE_n > 0 at every resolution (ci uphill ÔåÆ flux ciÔåÆcj)

  res   dWSE_n (m)    expected (m)  rel_err (%)  sign    
  ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
  12    0.470822      0.470822      0.000       Ô£ô
  13    0.210558      0.210558      0.000       Ô£ô
  14    0.094164      0.094164      0.000       Ô£ô
  15    0.042112      0.042112      0.000       Ô£ô
  16    0.018833      0.018833      0.000       Ô£ô
  17    0.008422      0.008422      0.000       Ô£ô
  18    0.003767      0.003767      0.000       Ô£ô
  19    0.001684      0.001684      0.000       Ô£ô
  20    0.000753      0.000753      0.000       Ô£ô

  Correct sign (dWSE_n > 0): 9 / 9
  Relative error: mean=0.000%  max=0.000%  (threshold: 2%)
Test Summary:               | Pass  Total  Time
GD2 ÔÇö d╠é projection formula  |   24     24  0.1s
ÔöÇÔöÇ GD3: n╠é_f projection formula ÔÇö tilted edges, multi-resolution ÔöÇÔöÇ
   Checks that dWSE_n recovers the face-normal WSE gradient component
   for both tilt directions at each A5 resolution 12ÔÇô20

  res    tilt(┬░)   dWSE_n (m)    expected (m)  rel_err(%)    sign    
  ÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇÔöÇ
  12     -25.0     0.426710      0.426710      0.000         Ô£ô
  12     25.0      0.426710      0.426710      0.000         Ô£ô
  14     -25.0     0.085342      0.085342      0.000         Ô£ô
  14     25.0      0.085342      0.085342      0.000         Ô£ô
  16     -25.0     0.017068      0.017068      0.000         Ô£ô
  16     25.0      0.017068      0.017068      0.000         Ô£ô
  18     -25.0     0.003414      0.003414      0.000         Ô£ô
  18     25.0      0.003414      0.003414      0.000         Ô£ô
  20     -25.0     0.000683      0.000683      0.000         Ô£ô
  20     25.0      0.000683      0.000683      0.000         Ô£ô
Test Summary:                 | Pass  Total  Time
GD3 ÔÇö n╠é_f projection formula  |   20     20  0.2s

All GD tests complete.
