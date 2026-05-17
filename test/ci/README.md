# QAtlas test architecture — design guardrails

The suite is organised around **three separated planes**.  Most CI bugs
historically came from blurring them, so the separation is enforced
structurally, not by convention.  Read this before changing anything
under `test/ci/`, `test/util/verify.jl`, `runtests.jl`, or `CI.yml`.

```
WHAT runs    →  test/ci/universe.jl     (single source of truth + guard)
HOW it splits → test/ci/plan_shards.jl  (timing-balanced; advisory only)
WHY it's right → test/util/verify.jl    (black-box cards → ci-evidence)
```

The golden rule: **WHAT is authoritative; HOW and WHY may never override
it.**  Timing/evidence data can be stale, missing, or wrong and the
suite must still run every test exactly once.

---

## WHAT — the test universe

- `test/ci/universe.jl` `ALL_DIRS` is the **only** place the test
  directory set is declared.  Do **not** re-list directories in
  `CI.yml`, the planner, or anywhere else.
- A test file **must** be named `test_*.jl` (the glob filter).  A file
  that does not match is silently ignored — name it correctly.
- Adding a new test directory **requires** adding it to `ALL_DIRS`.
  The completeness guard hard-fails CI (in every shard) if any on-disk
  `test_*.jl` directory is missing from `ALL_DIRS`, or if an `ALL_DIRS`
  entry is missing/empty.  This failure is **intended** — it makes a
  silently-skipped test structurally impossible.
- `test/test_aqua.jl` at the test/ root is the only sanctioned
  non-`ALL_DIRS` file (run once, in exactly one shard).

## HOW — sharding, timing, profiles

- Shard count `N` lives **only** in the `CI.yml` `plan` job.  Change
  parallelism there and nowhere else.
- `ci-timings` / `ci-evidence` are **orphan branches written only by
  `push:main`**; PR runs read them read-only.  Never write them from a
  PR (push contention + diff pollution).  Missing/stale → graceful
  round-robin / "pending" fallback, never a leak.
- Selection precedence: `QATLAS_TEST_FILES` > `QATLAS_TEST_SHARD` >
  `QATLAS_TEST_GROUP` > all.  `QATLAS_TEST_FILES` may only name files
  in the universe (the planner cannot smuggle in non-globbed files).
- Profiles (`QATLAS_TEST_PROFILE`): `fast` = PR merge gate (small N,
  loose tol, **no emit**); `full` = `push:main` (deeper, **emits**
  timing+evidence); `nightly` = cron (deepest).  The persisted numbers
  always reflect the heavier `full` run.
- **Dense ED is exponential.**  A spin-S chain is `(2S+1)^N`.  Hard-cap
  `N` per model at the feasible ceiling and pass explicit caps:
  `verify_profile_Ns(; fast=…, full=…, nightly=…)`.  Never let the
  profile knob push `N` past the dense-ED limit (spin-1 ⇒ N ≤ 8 =
  3^8; matches src `_MAX_ED_SITES_S1`).  *This was a real bug — 3^12
  dense eigen is ~TB of RAM.*
- `test/.ci-out/` is gitignored emit scratch — never commit it.

## WHY — black-box verification (the epistemic core)

- **`src/` holds closed-form analytical values only.**  ED / numerical
  diagonalisation belongs in **tests** as an independent cross-check,
  never as the src implementation.  (PR #359 was closed for putting
  dense-ED finite-T into src.)
- Tests are **black-box**: they know only (a) the model's physics
  (its written Hamiltonian) and (b) the public
  `fetch(model, quantity, bc; kw…)` API.  A test must **not** read
  `QAtlas._internal`s, nor reuse a src matrix-builder for its
  independent route — rebuild the physics from `test/util/generic_ed.jl`
  so a src-builder bug cannot also corrupt the cross-check.
- `verify(...)`: the **subject is always `fetch(...)` computed inside
  the helper**.  There is no argument to re-type the src formula —
  keep it that way; a card can then never be circular.
- The independent route must be a genuinely different computation path.
  `route ∈ {:ed_finite_size, :sum_rule, :delegation_invariant,
  :limiting_case, :literature_value, :second_closed_form}`.  Do not add
  a route that smuggles the src derivation back in.
- One `verify` call = one card on hub `Model/Quantity/BC`.  Confidence
  is **network-style**: number of *independent* routes × their
  precision.  Prefer several different routes on one hub over one.
- The card JSONL schema is the contract for the future doc/forum and
  MCP.  Keep it stable and minimal — scalars only (subject, a short
  convergence vector, errors, refs); never dump full spectra/matrices.
- A new value runs in its own PR (assert-only); its timing/evidence is
  recorded only on the post-merge `push:main`.  It is "pending
  verification" until then — never silently absent.

---

## Mechanics that bite (operational caveats)

- **Stacked PRs**: #361 → #362 → #363 → #364 each builds on the prior
  as a safe fallback.  Merge in dependency order; never out of order.
- **Branch protection** required checks are matched by job *name* and
  are `["build", "All tests passed"]`.  Keep the aggregate gate job
  named exactly `All tests passed` — renaming/removing it strands
  branch protection and permanently blocks merges.  Per-shard job
  names are *not* required checks (so they may change freely).
- `GITHUB_TOKEN` pushes do **not** trigger workflows (intentional
  anti-recursion).  Relied on for `ci-timings`/`ci-evidence` (data,
  must not trigger CI) and the reason auto-format was moved from a
  committing bot to a check (`FormatFix` → `FormatCheck`, PR #355).
- Dev-host mechanics: use the absolute juliaup path in non-interactive
  SSH; `gh` is at `~/.local/bin/gh`; nested ssh/heredoc mangles
  Unicode (β ⟨⟩ ≤) and `$` — edit via a scp'd file or base64, not
  inline heredocs.  The worktree test env needs
  `Pkg.develop(path=pwd()); Pkg.instantiate()`, which pollutes
  `test/Project.toml` with an absolute `[sources]` path — always
  `git checkout test/Project.toml` before committing.
