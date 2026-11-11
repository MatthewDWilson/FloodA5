<#
## Testing Sequence

Work through these in order — each gate confirms the previous layer before adding complexity. Paste the full terminal output for each stage back and I'll review before you proceed.

---

### Stage 1 — Unit tests (no mesh, no pya5 needed, ~10 seconds total)

Run these four existing test files first. They exercise the `EdgeList` constructor changes and confirm nothing is broken by the struct extension.

```powershell
#>

julia --project=. test\test_edge_geometry.jl
julia --project=. test\test_noc_correction.jl
julia --project=. test\test_sgs_unit.jl
julia --project=. test\test_inflow_point.jl

<#
```

**What to look for:** all tests pass with the same counts as before the change. Any `MethodError` on `EdgeList(...)` means a constructor site was missed. Any physics regression (sign wrong, mass not conserved) means the `dx_m`/`dy_m` values in a fixture are wrong.

---

### Stage 2 — New gradient direction test (no mesh, ~2 seconds)

```powershell
#>

julia --project=. test\test_gradient_direction.jl

<#
```

**What to look for:**

- GD1 table: `V̂_y` column should show **mixed signs** across res 12–20 (some positive, some negative). This confirms the diagnosed root cause is real on your machine's A5 geometry approximation.
- GD2 table: `sign` column should be **all ✓**, `rel_err (%)` should be **well below 2%** at every resolution.
- GD2c should report a CV well below 5%.

If GD2 shows any `✗ WRONG` in the sign column, stop and paste the output — that would indicate the formula has a sign error in the implementation.

---

### Stage 3 — Planar embankment re-run (two short runs, ~5–10 min each)

Use the same mesh and injection point as before, but now at res-16 first to confirm we haven't regressed what was working:

```powershell
# Baseline (unchanged — skip if you still have planar_result_std.h5 from before)
#>

julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh16_std.parquet `
    --flow-model standard `
    --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries `
    --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res16_baseline.h5 `
    --output-interval 300 `
    --gradient-correction off

# Corrected
julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh16_std.parquet `
    --flow-model standard `
    --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries `
    --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res16_corrected.h5 `
    --output-interval 300 `
    --gradient-correction on

<#
```

Then run the symmetry test:

```powershell
#>

julia --project=. test\test_planar_symmetry.jl `
    --baseline  test/planar_embankment/planar_res16_baseline.h5 `
    --corrected test/planar_embankment/planar_res16_corrected.h5 `
    --source-lat 51.0001 --source-lon -0.0434 `
    --sweep 10

<#
```

**What to look for:** the `|asym| gain` column should be **consistently positive** (corrected closer to zero than baseline) at most frames — same pattern as the previous Stage 3 result. The corrected asymmetry should stay on the **same side as zero** (not flip sign). Rough target: `|asymmetry|` corrected < 0.15 at the final frame.

---

### Stage 4 — Res-18 re-run (the failing case)

```powershell
# Baseline (skip if you still have planar_result_std.h5 from before)
#>


julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh18_std.parquet `
    --flow-model standard `
    --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries `
    --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res18_std_baseline.h5 `
    --output-interval 300 `
    --gradient-correction off

# Corrected
julia --project=. --threads auto FloodModel.jl `
    --meshload test/planar_embankment/planar_mesh18_std.parquet `
    --flow-model standard `
    --injection-point 51.0001,-0.0434,0.1 `
    --closed-boundaries `
    --sim-duration 72000 --dt-max 10 `
    --output test/planar_embankment/planar_res18_std_corrected.h5 `
    --output-interval 300 `
    --gradient-correction on


<#
```

```powershell
#>

julia --project=. test\test_planar_symmetry.jl `
    --baseline  test/planar_embankment/planar_res18_std_baseline.h5 `
    --corrected test/planar_embankment/planar_res18_std_corrected.h5 `
    --source-lat 51.0001 --source-lon -0.0434 `
    --sweep 10

<#
```

**What to look for:** this is the critical check. The corrected asymmetry should be **smaller in magnitude than baseline and not opposite in sign** — i.e., the plume should no longer flip from northward to southward. A corrected `|asymmetry|` below 0.3 at the final frame with the same sign as zero would be a clear pass. If the sign still flips (corrected asymmetry goes strongly negative), paste the output and the GD1 table from Stage 2 together — there's a further sign issue to diagnose.

---

### Stage 5 — Carlisle regression (confirm no stability regression)

Only proceed here if Stages 3 and 4 both look good.

```powershell
#>

julia --project=. --threads auto FloodModel.jl `
    --meshload test/carlisle/carlisle_mesh16_sgs.parquet `
    --rainfall 50 `
    --flow-model sgs `
    --sim-duration 36000 `
    --gradient-correction on `
    --output test/carlisle/carlisle_reg_corrected.h5 `
    --output-interval 600

<#
```

**What to look for from the run log:** `mb_err` stays at `0.0` throughout, `wet` cell count is stable, no oscillations in `max_depth`, `dt` doesn't collapse to near-zero. You don't need to paste the full log — just the final few checkpoint lines and the summary line at the end.

---

### What to paste back at each stage

| Stage | Paste |
|---|---|
| 1 | Full terminal output of all four test files |
| 2 | Full output of `test_gradient_direction.jl` |
| 3 | Full output of `test_planar_symmetry.jl` (res-16) |
| 4 | Full output of `test_planar_symmetry.jl` (res-18) + a screenshot of both Makie visualisations if convenient |
| 5 | Last ~20 lines of the Carlisle run log |

If anything fails at Stage 1 or 2, stop there — no point running the slower mesh-based stages until the unit tests are clean.
#>
