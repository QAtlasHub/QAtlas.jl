# Pre-PR checklist

Run these before opening a PR. CI enforces all of them; doing them locally
avoids a red PR.

## 1. Format (JuliaFormatter v2, pinned)

CI's FormatCheck runs `format(".", overwrite=false)` with **JuliaFormatter v2**.
Match it exactly — v1 formats differently:

```bash
julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add(name="JuliaFormatter", version="2"); \
          using JuliaFormatter; format(".")'
```

Style is `blue`, `margin=92` (see `.JuliaFormatter.toml`).

## 2. Regenerate the atlas + confirm the drift guard

Any change to a `*_registry.jl` (`@register`) **requires** regenerating the
atlas, or the inventory drift guard fails:

```bash
julia --project=. docs/atlas/generate.jl
julia --startup-file=no test/harness/atlas/test_inventory_drift.jl
julia --startup-file=no test/harness/atlas/test_doc_structure.jl
```

`generate.jl` is idempotent. Stage the substantive diffs (the new hub /
`ATLAS:DOCS` sections, `api.md`); see [documentation.md](documentation.md) for
the regen-noise hygiene.

## 3. Run the test suite

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Run the **multi-file** suite, not just your one file — a single-file green can
hide a util/method collision (see [verification.md](verification.md)).
`QATLAS_TEST_FILES="d/f.jl,…"` runs a subset; `QATLAS_TEST_FULL=1` enables the
heavy N=16 entanglement tests (nightly).

## 4. Build the docs locally

The Documenter build runs on **every** PR (not just `main`). Build it locally
and converge to green before pushing — broken `@ref`s and doctests fail the
build:

```bash
julia --project=docs docs/make.jl
```

If you removed / curated rendered docstrings, expect `@ref` breakage in a *wave*
(rendering one docstring exposes its own refs); iterate locally until clean.

## 5. Bump the version

Patch bump per PR in `Project.toml` (minor for larger additions). CI's
`check-version` requires the version to differ from `main`.

## 6. Open the PR

- Do **not** self-merge code PRs — the maintainer reviews and merges. Only
  trivial bot PRs may auto-merge.
- For a fix to an externally-reported issue, write `Refs #N`, not `Closes #N`
  (don't auto-close someone else's issue). Your own roadmap issues may be closed.

---

## Runtime gotchas

- **`QAtlas.fetch` vs `Base.fetch`** — always qualify `QAtlas.fetch(...)` in
  test code.
- **`fetch_kw` in `verify`** — pass parameters as a `NamedTuple`
  (`fetch_kw=(; beta=β)`). Passing them in the outer `verify(...; beta=β)` call
  silently swallows them.
- **OBC vs PBC gap** — in the ordered phase of the TFIM (`h ≪ J`) the lowest ED
  "gap" is the Z₂ tunnelling splitting (exponentially small in N), not the
  physical excitation gap.
- **Bond counting in small PBC lattices** — `Lattice2D.bonds(lat)` double-counts
  at `Lx=2` / `Ly=2`; the transfer-matrix and brute-force paths use the same
  convention so they agree, but be aware.
- **Entanglement-entropy RAM** — ED entanglement is `O(2^N)` memory; PR-CI runs
  `N ≤ 14`, `QATLAS_TEST_FULL=1` enables `N=16` (~24 GB), nightly only.
- **`git checkout main` with uncommitted work** — aborts; the changes carry over
  to the next branch. Commit or stash first.
- **PowerShell `Out-File`** — emits a BOM that breaks TOML/POSIX parsers; write
  via the editor or `utf8NoBOM`.
