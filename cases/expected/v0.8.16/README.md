# Expected values as of 0.8.16

The case metrics recorded up to 0.8.16, kept as the record when the SGS
source-order fix (`--source-order post`, the SGS default) re-recorded every
case with SGS runs. Standard-flow and SWE metrics in these files are
unchanged by it. See `CHANGELOG.md` and `docs/HYDRAULICS.md` §15.15.

To reproduce one, check out tag `v0.8.16` and run the case there
(`julia --project=. cases/run.jl run CASE`), or run the SGS jobs with
`--source-order pre`.
