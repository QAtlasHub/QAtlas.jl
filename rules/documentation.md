# Documentation

## Two layers

| Layer | Where | Written |
|-------|-------|---------|
| Framework reference | `docs/src/api.md` | by hand — Registry / Model / Quantity, scoped `@autodocs Pages=["core/…"]` |
| Per-model API | `docs/src/models/**.md` | **auto-generated** per-model `@autodocs` (see below) |
| Derivations | `docs/src/calc/*.md` | by hand — one topic per file, step-by-step |
| Atlas | `docs/src/atlas/**` | **auto-generated** by `docs/atlas/generate.jl` — never hand-edit |

`api.md` documents only the cross-cutting framework. `fetch(::Model, …)` only
exists as per-model `@register` + method pairs, so each model's `fetch`
docstrings live on **its own model page**, injected by the generator.

## Per-model `@autodocs` injection

`docs/atlas/generate.jl` scans `src/models` for each model's
`struct … <: AbstractQAtlasModel`, globs that model's source directory, and
writes an `@autodocs Pages=[…] Private=false` block between
`<!-- ATLAS:DOCS:START -->` / `<!-- ATLAS:DOCS:END -->` markers on the model
page — the model struct + every `fetch(::Model, …)` method + exported helpers,
in lock-step with `@register`. `generate.jl` is a **QAtlas-free static source
view** (it parses `*_registry.jl`, never loads the package), so it discovers
files by filesystem scan, not `methods(fetch)` introspection. Do not hand-edit
inside the `ATLAS:DOCS` (or `ATLAS:HUBS`) markers.

## The `@ref` web — the rule that bites

QAtlas docstrings are densely cross-linked with `` [`X`](@ref) ``. A symbol's
`@ref` resolves only if `X`'s docstring is rendered **somewhere** in the manual.
Consequences:

- **Do not add an `@autodocs`-everything dump.** It would render every internal
  `_*` helper. Render the framework on `api.md` (core files) and each model on
  its page; that covers the public surface without the dump.
- **`checkdocs=:none`** in `docs/make.jl`: the docs are curated per page, not a
  single index, so we do not require every export in the manual. Cross-references
  and doctests are still checked strictly.
- When you reference a symbol with `` [`X`](@ref) ``, make sure `X` is rendered:
  a public framework symbol → `api.md`; a model symbol → its model page. An
  **internal / non-exported** symbol (`_bc_size`, `pfaffian`, a `_model_*`
  helper, a Lattice2D type like `Honeycomb`) should be plain `` `X` ``, **not**
  an `@ref` — otherwise the build fails with an unresolved cross-reference.
- `size_threshold = 1_000_000` in `docs/make.jl` guards against an accidentally
  huge generated page. The per-page split keeps every page well under it.

## Atlas regeneration is mandatory after `@register`

```bash
julia --project=. docs/atlas/generate.jl
```

`generate.jl` is **idempotent**: re-running it twice is byte-stable. It also
re-touches per-model markdown; when committing, stage the *substantive* diffs
(your model page's `ATLAS:DOCS` / hub section, `api.md`) and let `git add docs/`
normalise CRLF — the atlas pages use plain links (no `@ref`), so unrelated
regen-reconciliation churn can be reverted safely to keep a PR focused.

## Derivation notes (`docs/src/calc/`)

One topic per file (`<math-object-slug>.md`, kebab-case). Section skeleton, in
order: `## Main result` (boxed result + 1-line validity), `## Setup`,
`## Calculation` (ending with a **Limiting-case checks** subsection that
evaluates the result at 2–3 canonical points against literature),
`## References`, `## Used by` (≥1 backlink). Bidirectional links are required:
the consumer page links to the derivation, the derivation lists its consumers.

**Depth standard — the calculation must be reconstructable.** Each line follows
from the previous by one *named* transformation (definition, algebraic identity,
change of variable, residue at a stated pole, …). Evaluate integrals (choose
the contour, list the poles, compute the residues); only genuinely external
theorems (Szegő, Cauchy, Stone–Weierstrass) may be cited, and then state the
theorem in one line before applying it. These phrases are **forbidden** (CI
greps for them, zero matches):

> "it can be shown", "we omit details", "standard calculation gives",
> "as in [textbook]", "one can verify", "it is easy to see", "it follows immediately".

If it takes one line, write the line. If it takes a page, write the page.

## Language

`docs/src/` and `rules/` are **English only** (public). Citations: `Author,
Journal Vol, page (Year), Eq. (n).` with a year on every entry (see
[citations.md](citations.md)).
