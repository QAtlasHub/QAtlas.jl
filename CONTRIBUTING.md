# Contributing to QAtlas.jl

Thank you for your interest in contributing to QAtlas! This guide describes how the project is organized and what we look for in contributions.

## What is QAtlas?

QAtlas is a **dictionary of rigorous results in quantum and statistical physics**. It stores analytically-known exact values in `src/` and cross-validates them against independent numerical / closed-form sources in `test/`. The core value proposition is that every stored result is checked against at least one *theoretically independent* derivation.

> **The checklist-grade implementation rules live in [`rules/`](rules/README.md).** Three non-negotiables:
> 1. **Every value is independently verified** — a `verify(...)` card reproduces it by a theoretically distinct route ([rules/verification.md](rules/verification.md)).
> 2. **Every value cites a precise DOI, checked against the paper** — download it with `doiget` and confirm the published value in the paper's own conventions before implementing ([rules/citations.md](rules/citations.md)).
> 3. **Regenerate the atlas + format before every PR** ([rules/pre-pr-checklist.md](rules/pre-pr-checklist.md)).

## Design Principles

### `src/` is a leaf — no lattice-package dependencies

The source code in `src/` does not depend on Lattice2D, QuasiCrystal, or any lattice construction package. It contains pure functions that map `(model, quantity, bc) → value`. Lattice-package dependencies live exclusively in `test/` via `[extras]`.

### Accumulate results first, refactor later

Adding a new rigorous result is always more valuable than perfecting the code structure. If a new result doesn't fit cleanly into the existing layout, **add it anyway and verify it** — the structure can be refactored later without losing the verified result.

### Physical correctness is paramount

A value in `src/` is only considered rigorous once it has been **independently verified** in `test/` via a `verify(...)` card (see [Verification cards](#verification-cards) below). Internal consistency checks (e.g. scaling relations) are necessary but not sufficient.

### Every numerical value traces to a derivation

For each new rigorous result, the accompanying `docs/src/calc/` derivation must be complete and step-by-step (see [rules/documentation.md](rules/documentation.md)). Forbidden phrases such as "it can be shown" / "we omit details" / "standard calculation gives" must not appear.

## Repository Structure

```text
src/
├── QAtlas.jl                          # Top-level module, includes everything
├── core/                              # Type hierarchy, registry, atlas plumbing
│   ├── type.jl                        #   AbstractQAtlasModel, OBC(N), PBC(N), Infinite, Quantity{S}
│   ├── quantities.jl                  #   Energy, FreeEnergy, MagnetizationX, SusceptibilityZZ, ...
│   ├── registry.jl                    #   @register macro + Implementation row schema
│   └── conversions.jl                 #   Generic <-> per-site Energy granularity dispatch
├── deprecate/                         # Pre-v0.13 Symbol-dispatch shims (removable at v1.0)
├── universalities/                    # Universality{C}(d) parametric type
│   ├── Universality.jl                #   base type + exponent table
│   ├── Ising2D.jl                     #   2D/3D Ising universality
│   ├── ONModel.jl                     #   XY / Heisenberg / O(N) σ
│   └── ...
└── models/
    │                                  # Layout: <class>/<Model>/<Model>.jl
    │                                  #         <class>/<Model>/<Model>_registry.jl
    ├── classical/
    │   ├── IsingSquare/
    │   │   ├── IsingSquare.jl         #   IsingSquare(; J, Lx, Ly) + fetch dispatches
    │   │   └── IsingSquare_registry.jl#   @register declarations
    │   └── IsingTriangular/
    │       ├── IsingTriangular.jl
    │       └── IsingTriangular_registry.jl
    └── quantum/
        ├── TFIM/                       #   Multi-file model: TFIM_thermal.jl, TFIM_dynamics.jl, ...
        ├── XXZ/
        ├── Hubbard1D/
        └── TightBinding1D/

test/
├── runtests.jl                         # Entry: respects QATLAS_TEST_PROFILE = fast|full|nightly
├── util/                               # Reusable test helpers
│   ├── verify.jl                       #   verify(...) card framework — READ THIS FIRST
│   ├── bloch.jl                        #   Generic Bloch Hamiltonian builder
│   ├── classical_partition.jl          #   Brute-force partition function
│   ├── tight_binding.jl                #   Real-space TB Hamiltonian
│   ├── spinhalf_ed.jl                  #   Dense spin-1/2 many-body ED
│   ├── sparse_ed.jl                    #   Sparse ED + KrylovKit Lanczos
│   └── extrapolate.jl                  #   1/N → ∞ extrapolation helpers
├── util_verify/                        # Self-tests of the verify framework
├── harness/atlas/                      # Atlas inventory drift guard, evidence card schema tests
│   └── test_inventory_drift.jl         #   Regen + diff guard for INVENTORY.jsonl
├── models/                             # Per-model verify-card test files
│   ├── classical/
│   └── quantum/
├── INVENTORY.jsonl                     # Frozen inventory of registered hubs (regenerated, drift-guarded)
└── verification/                       # Cross-model checks

docs/
├── atlas/
│   └── generate.jl                     # Regenerates docs/src/atlas/* from REGISTRY + INVENTORY
└── src/                                # Documenter source
    ├── calc/                           # Step-by-step derivations (Zettelkasten)
    ├── conventions.md                  # Sign / S / occupation conventions — REQUIRED READ
    └── atlas/                          # Auto-generated hub pages (do NOT hand-edit)

rules/                                  # Checklist-grade implementation rules — READ THESE
├── citations.md                        #   DOI policy + doiget literature cross-check
├── verification.md                     #   verify cards, independent routes
├── registry-conventions.md             #   @register, status axis, scheme=, CONVENTION header
├── documentation.md                    #   per-model @autodocs, generate.jl, @ref web, derivation depth
└── pre-pr-checklist.md                 #   format, atlas regen, version, docs build, gotchas
```

## The model API

Every model is a **concrete struct** with typed physical parameters; quantities and boundary conditions are likewise concrete structs (not `Symbol`s):

```julia
fetch(TFIM(; J=1.0, h=1.0), Energy(), OBC(24); beta=5.0)
fetch(XXZ1D(; J=1.0, Δ=0.5), LuttingerParameter(), Infinite())
fetch(IsingSquare(; J=1.0, Lx=4, Ly=4), PartitionFunction(); β=0.44)
```

Legacy Symbol calls still work through `src/deprecate/` but emit a deprecation log.

## How to Contribute

### Adding a new rigorous result

The minimum viable contribution is **four things** for a single `(model, quantity, bc)` triple:

1. **A `fetch(...)` dispatch** in `src/models/<class>/<Model>/<Model>.jl` (or a sibling file under the same directory).
2. **A `@register(...)` declaration** in `src/models/<class>/<Model>/<Model>_registry.jl`.
3. **A `CONVENTION` header** at the top of the source file (enforced by `test/lint/`; see `docs/src/conventions.md`).
4. **At least one `verify(...)` card** in `test/models/<class>/<Model>/...` that checks the value against an independent source.

#### 1. The `fetch` dispatch

```julia
"""
    fetch(m::MyModel, ::Energy{:per_site}, ::Infinite;
          beta::Real, J=m.J, kwargs...) -> Float64

Per-site internal energy at inverse temperature `β > 0`, ... (cite the
specific equation/theorem from the reference).

# References
- Author, *Journal* **Vol**, page (Year), Eq. (X.Y).
"""
function fetch(
    m::MyModel,
    ::Energy{:per_site},
    ::Infinite;
    beta::Real,
    J::Real=m.J,
    kwargs...,
)
    beta > 0 || throw(DomainError(beta, "MyModel Energy requires β > 0; got β = $beta."))
    return -J * tanh(beta * J)        # Author Year, Eq. (X.Y)
end
```

Conventions:

- Boundary condition is the third positional argument (`::OBC`, `::PBC`, `::Infinite`), **never** a kwarg.
- Physical parameters come from the model struct with kwarg overrides (`J::Real=m.J`).
- Validate inputs at the function boundary and throw `DomainError` with a precise message.
- Cite the source **in the docstring** to a specific equation/theorem — not just "Author (Year)".

#### 2. The registry declaration

`src/models/<class>/<Model>/<Model>_registry.jl`:

```julia
@register(
    MyModel,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/<class>/<Model>/test_<model>.jl",
    references=["Author Year"],
    notes="u(β,h=0) = -J tanh(βJ) per site; closed form, no integrator.",
)
```

| Field         | Allowed values                                                                                       |
| ------------- | ---------------------------------------------------------------------------------------------------- |
| `method`      | `:analytic`, `:bdg`, `:dense_ed`, `:sparse_ed`, `:transfer_matrix`, `:bethe_ansatz`, `:tba`, `:pfaffian`, `:not_implemented` |
| `reliability` | `:high` (closed form + literature-tested), `:medium` (ED only / cross-check), `:low` (heuristic), `:not_implemented` |
| `tested_in`   | Relative path to the test file that verifies this triple, or `nothing`                               |
| `references`  | `["Author Year", ...]` — keep short, point to the docstring for full details                         |
| `notes`       | One-line caller-facing caveat                                                                        |

The registry is consumed by `docs/atlas/generate.jl` to auto-generate `docs/src/atlas/hubs/<Model>_<Quantity>_<BC>.md` and the inventory drift guard.

A row also carries a `status` (the 4-value axis `:exact` / `:approx` / `:bound` / `:universal`, orthogonal to `reliability`). A `:approx` row (e.g. a high-temperature series) **requires** a `valid_domain` and usually an `error_order`, and is a second `scheme=…, canonical=false` definition of a hub whose canonical row stays exact — `fetch(m, q, bc; scheme=…)` selects it. See [rules/registry-conventions.md](rules/registry-conventions.md) for the `status`/`scheme=` rules and the `@eval`-method-loop caveat.

#### 3. The CONVENTION header

Every model source file must declare its sign / spin / occupation convention near the top:

```julia
# CONVENTION
#   Hamiltonian: -J Σ σ^z σ^z - h Σ σ^x   (FM convention, h ≥ 0)
#   Observable:  Spin-1/2 (S = σ/2)
#   Reference:   docs/src/conventions.md §Spin convention
```

The CI lint in `test/lint/` will reject merges that omit this. See [docs/src/conventions.md](docs/src/conventions.md) for the project-wide convention policy and which Hamiltonian sign / observable normalization applies to which model family.

### Verification cards

**Every new `(model, quantity, bc)` triple requires at least one `verify(...)` card.** A card black-box-checks `fetch(...)` against an `independent` value obtained by a route theoretically distinct from the implementation.

The signature (frozen — see `test/util/verify.jl`):

```julia
verify(
    model,
    quantity,
    bc;
    route::Symbol,                      # one of the allowed routes below
    independent,                        # scalar OR convergence vector
    agree_within::Real,                 # absolute tolerance
    refs::AbstractVector{<:AbstractString},
    fetch_kw::NamedTuple = (;),         # passed into the black-box fetch
    reliability::Symbol = :high,
    at = nothing,                       # for convergence vectors: the x-axis values
    expected_fail::Bool = false,        # @test_broken (bug-surfacer)
    subject_extract = nothing,          # project a NamedTuple/Vector fetch result to one Float64
) -> subject
```

Allowed routes (see `_VERIFY_ROUTES` in `test/util/verify.jl`):

| Route                          | When to use                                                                       |
| ------------------------------ | --------------------------------------------------------------------------------- |
| `:ed_finite_size`              | Exact diagonalisation at small `N`, optionally extrapolated to `N → ∞`            |
| `:second_closed_form`          | A different closed-form derivation of the same quantity                           |
| `:limiting_case`               | Known value at a special point (`T = 0`, `μ = 0`, `Δ = 0`, β → 0/∞ limit, ...)    |
| `:sum_rule`                    | Independent analytic identity (e.g. f-sum rule)                                   |
| `:delegation_invariant`        | Model X at parameter p is exactly model Y (e.g. `XXZ1D(Δ=0)` ≡ free fermion)      |
| `:literature_value`            | Published numeric cross-check (DMRG / MC table)                                   |
| `:lieb_square_ice`             | Lieb 1967 ice-point closed form (SixVertex)                                       |
| `:lieb_ferroelectric`          | Lieb 1967 frozen-FE closed form (SixVertex)                                       |
| `:single_root_specialisation`  | Specialised single-root closed form (e.g. SU(2)_k Verlinde S₀₀)                   |
| `:multi_root_product`          | Multi-positive-root product form (e.g. SU(N≥3)_k Verlinde S₀₀)                    |

Worked example — `IsingChain1D` finite-T entropy `s(β,h=0) = log(2 cosh βJ) − βJ tanh βJ`:

```julia
@testset "IsingChain1D — ThermalEntropy closed form" begin
    for (J, β) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 1.0), (2.0, 0.7))
        verify(
            IsingChain1D(; J=J),
            ThermalEntropy(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=log(2 * cosh(β * J)) - β * J * tanh(β * J),
            agree_within=1e-12,
            refs=["Ising 1925: s(β,h=0) = log(2 cosh βJ) − βJ tanh(βJ) per site"],
        )
    end
end
```

Critical properties:

- The `subject` is **fetched inside** `verify` — you cannot pre-compute it on the caller side and pass it in. This prevents tautological cards.
- The `independent` value must be derivable **without** running the same code path as `fetch`. Re-deriving the same integral with the same `quadgk` call is *not* independent.
- Cite the source of the `independent` value in `refs` precisely (paper + equation / page).
- For a **literature** `independent` value, get it from the paper itself: download with `doiget`, read the published number in the paper's own conventions, and anchor one clean coefficient to the code. A self-derivation cannot catch a convention error (spin normalisation, sign, per-site vs per-bond). See [rules/citations.md](rules/citations.md).
- Tolerance `agree_within` is **absolute**. For relative tolerance against a value `x`, pass `agree_within = abs(x) * rtol`.

Limit-only cards (when there is no second closed form) use `route=:limiting_case` and cite the limit:

```julia
# β → 0⁺ limit of free-fermion free energy: ω → -log 2 / β per site
verify(
    TightBinding1D(),
    FreeEnergy(),
    Infinite();
    route=:limiting_case,
    fetch_kw=(; beta=1e-3),
    independent=-log(2) / 1e-3,
    agree_within=abs(log(2) / 1e-3) * 2e-3,
    refs=["Mahan §1.3: free-fermion β → 0⁺ limit ω → -T log 2"],
)
```

Tests that probe **exception shape** (e.g. `DomainError` on `β ≤ 0`) and tests that probe **identities between several fetched values** (e.g. Gibbs `s = β(u − f)`) do not fit the per-quantity card schema and stay as plain `@test_throws` / `@test`. Use them sparingly and as a complement to — not a replacement for — verify cards.

### Atlas regeneration

After adding or changing any `@register`, regenerate the atlas:

```bash
julia --project=docs docs/atlas/generate.jl
```

This rewrites every auto-generated atlas surface, all derived from the fixed substrate (registry + INVENTORY + R1 assurance + ED_INFEASIBLE_MODELS + calc/*.md filenames):

| Auto-generated page | Content | Derived from |
|---|---|---|
| `docs/src/atlas/index.md` | top atlas + risk-linter + per-model breakdown table | registry + INVENTORY |
| `docs/src/atlas/ModelList.md` | top searchable catalog, one row per model, columns: Universality, #K, methods, assurance distribution, ED-feasibility, regimes | substrate-derived |
| `docs/src/atlas/models/<Model>.md` × 58 | per-model `Quantity × BC` matrix, **Convention** block (from `# CONVENTION` header), **Derivation notes** (matched calc/*.md), aggregated methods + refs | registry + INVENTORY + src file comment + calc filenames |
| `docs/src/atlas/quantities/<Quantity>.md` × 51 | inverse `Model × BC` matrix, methods aggregation, universality coverage, top references | registry + INVENTORY |
| `docs/src/atlas/hubs/<Model>_<Quantity>_<BC>.md` × 263 | per-hub card with `src` claim, corroboration cards table, reconstructed `verify(...)` call, **Derivation note** link, three-way back-links (Model, Quantity, Atlas) | registry + INVENTORY |
| `docs/src/atlas/by/{model,quantity,bc,level,mechanism,regime}.md` | 1D facet aggregators | INVENTORY |
| `docs/src/atlas/Bibliography.md` | all citations, deduplicated, with hub backlinks sorted by hub-count | registry refs |
| `docs/src/atlas/CalcIndex.md` | inverse view: every `docs/src/calc/*.md` ↔ matched models | calc filenames |
| `docs/src/atlas/Audit.md` | doc-health gap surface: missing CONVENTION headers, missing quantity Definitions, orphan calc notes, models with 0 hubs, INVENTORY card hubs without registry claim (split universality-class vs real) | substrate-derived |
| `docs/src/atlas/Methods.md` | solution-mechanism facet — every registered route grouped by method, with model × quantity counts | registry mechanism |

The `test/INVENTORY.jsonl` drift guard then enforces that the regenerated inventory matches the committed one:

Adding a new `@register` entry therefore automatically:

1. creates a new per-hub card under `hubs/`,
2. adds the `(Quantity, BC)` cell to that model's `Quantity × BC` matrix in `models/<Model>.md` (empty cells are gap visualisation — where physics could be added next),
3. adds the `(Model, BC)` cell to that quantity's `Model × BC` matrix in `quantities/<Quantity>.md`,
4. bumps the `#K` count and assurance distribution in `ModelList.md`.

Two CI guards enforce that the auto-generated structure stays consistent with the substrate: `test/harness/atlas/test_inventory_drift.jl` (registry/INVENTORY drift) and `test/harness/atlas/test_doc_structure.jl` (per-model / per-quantity / ModelList structural completeness + per-hub back-link presence):

```bash
julia --startup-file=no test/harness/atlas/test_inventory_drift.jl
```

A PR that adds registrations without regenerating the atlas will fail this drift guard.

### Adding a universality class

```julia
fetch(Universality(:MyClass), CriticalExponents(); d=2)
```

- Use `Rational{Int}` for exact values (e.g. `β = 1//8`).
- For numerical estimates, include `_err` fields (e.g. `β = 0.32642, β_err = 1e-5`).
- Verify scaling relations: `α + 2β + γ = 2` (Rushbrooke), `γ = β(δ−1)` (Widom).

### Adding a tight-binding lattice

The generic Bloch builder in `test/util/bloch.jl` works for any Lattice2D topology:

```julia
λ_bloch = bloch_tb_spectrum(MyTopology, Lx, Ly, t)
H = build_tight_binding(lat, t)
λ_real = sort(eigvals(Symmetric(H)))
@test λ_bloch ≈ λ_real atol=1e-10
```

**Topology-name collisions with Lattice2D.** `Honeycomb`, `Kagome`, `Lieb`, `Triangular` exist in both packages and are **not exported** by QAtlas. Qualify as `QAtlas.Honeycomb()` etc. The alias `Graphene = Honeycomb` *is* exported (no collision).

## Documentation

Writing a new `docs/src/calc/` note? First read [rules/documentation.md](rules/documentation.md). The depth standard is enforced by:

- grep check: zero matches for `it can be shown`, `we omit`, `standard calculation`, `one can verify`, `it is easy to see`, `it follows immediately`.
- Structure: `## Main result`, `## Setup`, `## Calculation`, `## References`, `## Used by` — in that order.
- Documenter build: `julia --project=docs -e 'include("docs/make.jl")'` completes with zero `Error:` lines.

Exemplars: [docs/src/calc/jw-tfim-bdg.md](docs/src/calc/jw-tfim-bdg.md), [docs/src/calc/bethe-ansatz-heisenberg-e0.md](docs/src/calc/bethe-ansatz-heisenberg-e0.md).

Public docs (`docs/src/`), `rules/`, and `CONTRIBUTING.md` are **English only**.

`docs/src/atlas/*` is auto-generated — **do not hand-edit**. To change the rendering, edit `docs/atlas/generate.jl` and regenerate.

### Docs architecture (api.md + per-model `@autodocs`)

`docs/src/api.md` is a hand-written **framework reference** (Registry / Model / Quantity, scoped to the `core/` source files). It is *not* a full symbol dump. Each model's `fetch(::Model, …)` docstrings are rendered on **that model's own page**, via an `@autodocs` block that `docs/atlas/generate.jl` injects between `<!-- ATLAS:DOCS:START/END -->` markers (scoped to the model's source directory). So `fetch` documentation stays in lock-step with `@register`, on the model where it lives.

Because QAtlas docstrings are densely cross-linked with `` [`X`](@ref) ``, `make.jl` sets `checkdocs=:none` (curated per page, not one index) but still checks cross-references and doctests strictly. When you write `` [`X`](@ref) ``, ensure `X` is rendered somewhere — a framework symbol on `api.md`, a model symbol on its page; an **internal / non-exported** symbol must be plain `` `X` ``, not an `@ref`. See [rules/documentation.md](rules/documentation.md).

## CI gates

A PR must pass:

1. **Format check** — `JuliaFormatter.format(path)` is idempotent on every changed `.jl` file.
2. **Convention lint** — every modified model file has its `CONVENTION` header (`test/lint/`).
3. **Test suite** — `Pkg.test()`, possibly under different `QATLAS_TEST_PROFILE` profiles.
4. **Inventory drift guard** — `test/harness/atlas/test_inventory_drift.jl` passes (regenerate the atlas if it doesn't).
5. **Aqua** (`test/test_aqua.jl`) — no stale deps, no piracy, no ambiguities.
6. **Documenter build** — `julia --project=docs docs/make.jl` succeeds (run on every PR, not just on `main`).

## Things to watch out for

- **`QAtlas.fetch` vs `Base.fetch`**: always qualify as `QAtlas.fetch(...)` in test code.
- **Bond counting in small PBC systems**: Lattice2D's `bonds(lat)` double-counts when `Lx = 2` or `Ly = 2` with periodic boundaries. Both the transfer-matrix and brute-force paths use the same convention, so they agree — but be aware.
- **OBC vs PBC for gap analysis**: in the ordered phase of the TFIM (`h ≪ J`), the lowest ED "gap" is the Z₂ tunneling splitting (exponentially small in N), not the physical excitation gap.
- **ForwardDiff compatibility**: relax `Float64` → `Real` and avoid LAPACK-only operations (`eigvals(Symmetric(T))`) — use `tr(T^n)` with dual numbers.
- **Entanglement tests and RAM**: ED entanglement scales as `O(2^N)` memory at `N ≥ 14`. PR-CI runs `N ≤ 14`; `QATLAS_TEST_FULL=1` enables `N = 16` (~24 GB peak), run nightly by [`NightlyFullTests.yml`](.github/workflows/NightlyFullTests.yml).
- **`fetch_kw` in `verify`**: pass parameters as a `NamedTuple` (`fetch_kw=(; beta=β)`), not as positional. A common bug is passing them in the outer `verify(...; beta=β)` call where they're silently swallowed.

## Before submitting a PR

1. Run the full test suite locally:

   ```bash
   julia --project=. -e 'using Pkg; Pkg.test(; julia_args=`-t auto --heap-size-hint=96G`)'
   ```

2. If you added registrations, regenerate the atlas and confirm both guards:

   ```bash
   julia --project=docs docs/atlas/generate.jl
   julia --startup-file=no test/harness/atlas/test_inventory_drift.jl
   julia --startup-file=no test/harness/atlas/test_doc_structure.jl
   ```

3. Bump the version in `Project.toml` (patch bump per PR; minor for larger additions).

4. If you added a `docs/src/calc/` derivation, build the docs:

   ```bash
   julia --project=docs -e 'include("docs/make.jl")'
   ```

5. If you find a discrepancy between QAtlas values and the literature, open an issue with the full reference (journal, volume, page, equation number).
