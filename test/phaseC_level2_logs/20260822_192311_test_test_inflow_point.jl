[ Info: A5Grid: CUDA GPU detected ÔÇö PIP sampling will run on GPU
Ôöî Info: Starting FloodA5 model...
Ôöö   Dates.now() = 2026-08-22T19:24:59.118
Ôöî Info: FloodA5 model finished.
Ôöö   Dates.now() = 2026-08-22T19:25:00.466
Test Summary:                            | Pass  Total  Time
T-IP1: _interp_hydrograph at knot points |    3      3  0.7s
Test Summary:                             | Pass  Total  Time
T-IP2: linear interpolation between knots |    3      3  0.0s
Test Summary:                               | Pass  Total  Time
T-IP3: flat extrapolation before first knot |    2      2  0.0s
Test Summary:                             | Pass  Total  Time
T-IP4: flat extrapolation after last knot |    2      2  0.0s
Test Summary:                            | Pass  Total  Time
T-IP5: apply_source! adds correct volume |    2      2  0.1s
Test Summary:                                    | Pass  Total  Time
T-IP6: cumulative_volume trapezoidal integration |    4      4  0.1s
[ Info: BDY parse: file header: 'Test BDY header'
Test Summary:                          | Pass  Total  Time
T-IP7: LisfloodBDYReader hoursÔåÆseconds |    6      6  0.6s
[ Info: BDY parse: file header: 'Test BDY header'
Test Summary:                       | Pass  Total  Time
T-IP8: LisfloodBDYReader two series |    4      4  0.0s
[ Info: BDY parse: file header: 'Test BDY header'
Test Summary:                    | Pass  Total  Time
T-IP9: parse_bci_file QVAR entry |    8      8  0.3s
Test Summary:                     | Pass  Total  Time
T-IP10: parse_bci_file QFIX entry |    4      4  0.2s
Test Summary:                                      | Pass  Total  Time
T-IP11: BCI N/E/S/W entries log warning, not error |    2      2  0.0s
Test Summary:                                     | Pass  Total  Time
T-IP12: two InflowPoints on same cell sum volumes |    1      1  0.0s
[ Info: test_inflow_point.jl: all tests complete
