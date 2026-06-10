# QAtlas implementation rules

Focused, checklist-grade rules for contributing physics to QAtlas.
`CONTRIBUTING.md` is the narrative guide; this directory is the detail you
actually have to get right.

## Three non-negotiables

1. **Every stored value is independently verified.** A number in `src/` is not
   "rigorous" until a `verify(...)` card in `test/` reproduces it by a
   *theoretically distinct* route — exact diagonalisation, a second closed
   form, a limiting case, or a cross-model exact match. Internal consistency
   (scaling relations, autodiff identities) is necessary but **not** sufficient.
   → [verification.md](verification.md)

2. **Every value cites a precise DOI, checked against the paper itself.** Not
   "Author (Year)" — the exact DOI plus the equation/table number. Download the
   source with `doiget` and confirm the published value *in the paper's own
   conventions* (spin normalisation, sign, per-site vs per-bond, …) **before**
   you implement. A self-derivation can be internally consistent and still be
   wrong by a convention factor; only the paper settles it.
   → [citations.md](citations.md)

3. **Regenerate + format before every PR.** Run JuliaFormatter (v2, pinned),
   regenerate the atlas, and confirm the inventory drift guard. Adding an
   `@register` without regenerating the atlas fails CI.
   → [pre-pr-checklist.md](pre-pr-checklist.md)

## The rules

| File | Covers |
|------|--------|
| [citations.md](citations.md) | DOI policy + the `doiget` literature cross-check |
| [verification.md](verification.md) | verify cards, independent routes, no circular checks |
| [registry-conventions.md](registry-conventions.md) | `@register`, the 4-`status` axis, `scheme=`, `# CONVENTION` headers |
| [documentation.md](documentation.md) | per-model `@autodocs`, `generate.jl`, the `@ref` web, derivation-depth standard |
| [pre-pr-checklist.md](pre-pr-checklist.md) | format, atlas/inventory regen, version bump, docs build, runtime gotchas |

## Adding a result (quick reference)

For a single `(model, quantity, bc)` triple:

1. A `fetch(...)` dispatch + docstring citing the **DOI and equation**.
2. An `@register(...)` row (`method`, `status`, `reliability`, `references`, `tested_in`).
3. A `# CONVENTION` header on the source file.
4. At least one `verify(...)` card against an **independent** source.

Then: regenerate the atlas, format, bump the version, build the docs locally.
Each step is detailed in the rule files above.
