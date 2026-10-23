# FloodA5 — Versioning Policy

FloodA5 did not have formal version numbers before this document was
written. This is the policy adopted going forward, and the reasoning
behind it, recorded once here rather than re-litigated at each release.

## Scheme: Semantic Versioning, starting at 0.1.0

`MAJOR.MINOR.PATCH`, per [semver.org](https://semver.org):

- **PATCH** — any change that leaves **default results** unchanged: a bug
  fix that does not alter default output, a new *opt-in* option or
  experimental mode, diagnostics, tests, documentation, or a performance
  change with bit-identical results. (Since 0.2.1; before that PATCH meant
  a pure bug fix only. Under semver's `0.x` rule this is a project choice,
  made so that the steady stream of opt-in research options does not force
  a MINOR release each time.)
- **MINOR** — a change to default behaviour or results, a change to the
  output format, or the removal of an option. Under semver's own rule for
  the `0.x` series (see below), a MINOR bump *may* include breaking changes
  while the project is pre-1.0 — this is expected and fine, not a deviation
  from the standard. MINOR versions get a GitHub release.
- **MAJOR (1.0.0)** — reserved for the point at which the model is
  considered validated against real-world data and the CLI/output
  interfaces are considered stable enough that external users should
  expect breaking changes only at major version boundaries. **Not yet
  reached.** See the "Why 0.x, not 1.0" section below.

The current version is recorded in `VERSION` (repository root, plain text,
single line) and mirrored in `FloodModel.jl` as `const FLOODA5_VERSION`.
Keep both in sync — `VERSION` is what a build/release script would read;
the Julia constant is what gets printed by `--version` and the startup
banner. `julia FloodModel.jl --version` prints the running version.

## Why 0.x, not 1.0

Semver's own specification is explicit that the `0.x` series means
"anything may change at any time" — appropriate for software still under
active development, without a stability contract yet. That describes
FloodA5's actual state accurately:

- The directional-bias corrections are the default since 0.2.0 and have been
  checked on real terrain with open boundaries, but the standard solver still
  carries grid-scale noise in deep channels and its results are sensitive to
  the timestep (see `docs/HYDRAULICS.md` §15) — the physics is not yet pinned
  down.
- Default behaviour and CLI flags have changed multiple times across the
  project's development (sill definitions, CFL formula, boundary-condition
  defaults, flux-limiter behaviour) — the kind of change semver expects to
  see reflected honestly in a sub-1.0 version number, not hidden behind a
  1.0 that implies a stability promise the project isn't ready to make.
- No formal validation against observed flood extent data (e.g. an NSE/KGE
  skill-score comparison) has been completed yet — see the roadmap in
  `README.md`.

1.0.0 is the right target once those are addressed, not before.

## Why start at 0.1.0 specifically, and not something reflecting the amount of work already done

It was tempting to pick a "bigger" starting number — the model already has
a substantial, tested feature set (two solvers, sub-grid sampling,
dynamic inflows, boundary conditions, multiple directional-bias
correction candidates) built up over many development sessions and
several complete feature branches. But a version number's job is to
signal *compatibility going forward from a starting point*, not to
encode how much work already happened before that point existed. 0.1.0
as a first tag doesn't understate the project; it just means "this is
where formal versioning begins," with normal room to grow via MINOR bumps
from here.

## Why the project's internal bug-numbering (`Bug 61`, etc.) is not the version number

The project's development history tracks a running, numbered list of
bugs found and fixed (recorded internally, referenced from `HYDRAULICS.md`
and elsewhere in general terms). It was suggested that this count could
become the PATCH number directly — e.g. `0.1.60` after the 60th fix. This
is **not recommended**, for a concrete reason: PATCH, under semver, means
*only* a bug fix with no behavioural change. Many of the numbered "bugs"
in this project's history were substantial behavioural or interface
changes in their own right — adding the Froude and volume limiters,
changing the default sill definition, adding whole new CLI flags as part
of a fix. Encoding that count as PATCH would misrepresent those changes
to anyone relying on semver's actual meaning (e.g. tooling that assumes
PATCH-only upgrades are safe to apply automatically). The bug log remains
exactly where it already lives — the project's internal development
history — as a detailed engineering record; the public version number
tracks release-facing compatibility instead, using the ordinary semver
rules above.

## No retrofitted version numbers for prior development

Earlier commits, branches, and sessions are not being assigned version
numbers after the fact. Two reasons: first, there's no principled way to
decide exactly where a boundary between two fictional version numbers
should fall in history that was never tagged at the time — any answer
would be invented, not recovered. Second, the project already has a
complete, more useful historical record for this purpose — the internal
development history and `PROJECT_STATE.md`'s cumulative bug log — and
duplicating that as a sequence of fake version tags would create the
appearance of precision that isn't real. `CHANGELOG.md` instead opens
with a short, honest prose summary of the pre-0.1.0 development
trajectory (the major branches/phases, not individual bug fixes), clearly
separated from the versioned entries that follow it.

## Practical guidance for the next release

When you're ready to cut the next version:

1. Decide PATCH vs. MINOR using the rules above — the test is whether
   default results change. When in doubt, treat it as MINOR (semver's
   guidance: "if in doubt, bump the more significant number").
2. Every merged PR bumps the version (usually PATCH). Update `VERSION` and
   `FLOODA5_VERSION` together in the PR's **final** commit, not when the
   branch starts, so parallel branches do not conflict on the version.
3. Give the version its own `CHANGELOG.md` entry (`## [0.2.1] — date`);
   `[Unreleased]` collects work on branches not yet merged.
4. Tag the merge commit (`git tag -a v0.2.1 …`) so every version string is
   traceable to its code. Publish a GitHub release (notes from the
   CHANGELOG entry) for MINOR versions only.
