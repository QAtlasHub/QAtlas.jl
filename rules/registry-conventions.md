# Registry and conventions

## `@register`

Every `(model, quantity, bc)` triple is declared in
`src/models/<class>/<Model>/<Model>_registry.jl`:

```julia
@register(
    MyModel, Energy{:per_site}, Infinite,
    method=:analytic,
    status=:exact,
    reliability=:high,
    references=["AuthorYear"],
    tested_in="test/models/<class>/<Model>/test_<model>.jl",
    notes="One-line caller-facing caveat.",
)
```

| Field | Values |
|-------|--------|
| `method` | `:analytic`, `:bdg`, `:dense_ed`, `:sparse_ed`, `:transfer_matrix`, `:bethe_ansatz`, `:tba`, `:pfaffian`, `:htse`, … |
| `status` | `:exact`, `:bound`, `:approx`, `:universal` (the 4-value axis; see below) |
| `reliability` | `:high` (closed form + literature-tested), `:medium` (ED / cross-check), `:low` |
| `references` | bibkeys in `docs/references.bib` (see [citations.md](citations.md)) |
| `tested_in` | path to the verify-card test file |

## The 4-value `status` axis

`status` is **orthogonal** to `reliability` — it is the *kind* of claim:

- `:exact` — an analytic closed form, verified as an equality. **Forbids**
  `valid_domain` / `error_order`.
- `:approx` — a domain-limited approximation (a high-temperature series, a
  large-N proxy). **Requires** `references` *and* a `valid_domain`; usually also
  an `error_order`.
- `:bound` — a model-independent or one-sided bound; requires a `direction`.
- `:universal` — derived by construction from a `Universality{C}` / `Bound{D}`
  class, not a per-model claim.

## Multiple definitions of one hub (`scheme=`)

A single `(model, quantity, bc)` hub can carry several definitions, selected by
`scheme=` at fetch time. The default `scheme=:canonical` stays the exact row; an
approximation is a second `canonical=false` row:

```julia
@register(
    AKLT1D, SpecificHeat, Infinite,
    scheme=:htse, method=:htse, status=:approx,
    valid_domain="βJ ≲ 0.4 (high temperature)",
    error_order="O((βJ)⁵)",
    canonical=false,
    references=["Lohmann2014", "AKLT1988"],
    tested_in="test/models/quantum/misc/test_aklt_htse.jl",
    notes="…",
)
```

The fetch routes on the keyword: `fetch(m, q, Infinite())` reproduces the exact
canonical row; `fetch(m, q, Infinite(); scheme=:htse, beta=…)` routes to the
approximation. Mirror the existing TFIM `scheme=:high_T` pattern. Use
`definitions(model, quantity, bc)` to list a hub's schemes.

### Watch out: `@eval` method loops

Some models define their thermal `fetch` methods through an `@eval` loop over a
`_MODEL_THERMAL_METHODS` tuple (e.g. TFIM). A literal `::SpecificHeat` grep will
**miss** these. Before adding a `fetch(Model, Quantity, …)` method, check for an
`@eval` / `_*_METHODS` loop in the model's files — adding a duplicate triggers a
method-overwrite warning and a registry `canon` clash.

## Coherence: `errors = 0` is the hard invariant

`src/core/coherence.jl` (C1–C9) structurally self-reports holes.
`coherence_report()` must hold `:error` findings to **zero**. `:gap` findings
(self-reported coverage holes, e.g. orphan `:Potts3` / `:Potts4` classes) **need
not be empty** — they are honest "not done yet" markers, not failures. Coherence
is *structural*, not physical; the physics is checked by the verify cards.

## The `# CONVENTION` header

Every model source file declares its sign / spin / occupation convention near
the top — the CI lint rejects merges that omit it:

```julia
# CONVENTION
#   Hamiltonian: -J Σ σ^z σ^z - h Σ σ^x   (FM convention, h ≥ 0)
#   Observable:  Spin-1/2 (S = σ/2)
#   Reference:   docs/src/conventions.md §Spin convention
```

## `fetch` dispatch conventions

- Boundary condition is the **third positional argument** (`::OBC` / `::PBC` /
  `::Infinite`), never a kwarg.
- Physical parameters come from the struct with kwarg overrides (`J::Real=m.J`).
- Validate at the boundary and `throw(DomainError(x, "precise message"))`.
- Keep `src/` a leaf: **no** Lattice2D / QuasiCrystal dependency in `src/`
  (those live in `test/` via `[extras]`).
- `Honeycomb`/`Kagome`/`Lieb`/`Triangular` collide with Lattice2D and are **not
  exported** — qualify as `QAtlas.Honeycomb()`.
- For ForwardDiff-differentiated quantities, type parameters as `Real` (not
  `Float64`) and avoid LAPACK-only ops (`eigvals(Symmetric(·))`) — use forms
  that accept dual numbers.
