# Expected values before the SGS source-order fix

The case metrics recorded up to 0.8.x, kept as the record when the SGS
source-order fix (`--source-order post`, the SGS default) re-recorded every
case with SGS runs. Standard-flow and SWE metrics in these files are
unchanged by it. See `CHANGELOG.md` and `docs/HYDRAULICS.md` §15.15.

To reproduce one, check out the tag of that version and run the case there
(`julia --project=. cases/run.jl run CASE`), or run the SGS jobs with
`--source-order pre`.
