# Verification

QAtlas is a reference library: every stored value must be
**demonstrably correct**. The test suite is organised into five layers,
each with a distinct purpose and level of independence.

---

## Layer 1 — Unit tests (`test/core/`)

Verifies that the core machinery — registry, alias routing, type
hierarchy, energy granularity, Pfaffian implementation — behaves
correctly independently of any physical model.

| File | What is checked |
|------|-----------------|
| `test_registry.jl` | `@register` / lookup roundtrip; `implementation_status` |
| `test_alias.jl` | `OBC → PBC` and `PBC → Infinite` alias chains |
| `test_type.jl` | `Quantity` dispatch, `BoundaryCondition` hierarchy |
| `test_energy_granularity.jl` | `Energy{:total}` ↔ `Energy{:per_site}` conversion |
| `test_pfaffian.jl` | Pfaffian correctness for small matrices |
| `test_verify_harness.jl` | `verify()` card execution mechanics |

---

## Layer 2 — Identity and symmetry tests (`test/identities/`)

Checks model-specific mathematical identities using **only QAtlas
itself** — no external computation. Catches formula errors that
survive unit testing because the physics constraint is violated.

| File | What is checked |
|------|-----------------|
| `test_identities_TFIM.jl` | Scaling relations at criticality; $E_0$ BC consistency |
| `test_identities_TFIM_pbc.jl` | PBC parity-sector consistency |
| `test_identities_Heisenberg1D.jl` | SU(2) symmetry; $e_0$ bounds; dimer exact values |
| `test_identities_XXZ1D.jl` | Anisotropy limits $\Delta \in \{-1,0,1\}$; Luttinger $K, u$ |
| `test_identities_IsingSquare.jl` | Yang formula limits; $T_c$ self-consistency |
| `test_identities_KitaevHoneycomb.jl` | Flux-free sector ground state; gauge structure |
| `test_identities_S1Heisenberg1D.jl` | Haldane gap; AKLT point |
| `test_cross_bc_scaling.jl` | OBC → PBC → Infinite consistency across system size |
| `test_TFIM_limits_cross_model.jl` | $J=0$ / $h=0$ limits agree with IsingChain1D and free fermions |
| `test_TFIM_dynamic_symmetries.jl` | Time-reversal and parity in dynamic correlators |
| `test_property_invariants.jl` | Monotonicity, positivity, and bound constraints |

**What this layer cannot do**: It cannot detect systematic errors
where the source formula itself is wrong, because there is no
independent computation to compare against.

---

## Layer 3 — Physics cross-verification (`test/verification/`)

Cross-checks `src/` analytical formulas against **independent ED** via
`Lattice2D`-built real-space Hamiltonians, or against independent
integration (Yang-Yang). The two paths must agree to numerical precision.

| Sub-directory | Coverage |
|---------------|----------|
| `tfim_ising/` | BdG gap closure; FDT sanity and numerics; AD thermodynamics vs transfer-matrix |
| `heisenberg_xxz/` | Luttinger parameters (Bethe vs ED); Yang-Yang integral |
| `universality/` | Universality ↔ model cross-checks (8 connections); thermodynamic identities |

See [Cross-Verification Table](cross-checks.md) for all 8 universality ↔ model
connections.

---

## Layer 4 — Atlas harness (`test/harness/atlas/`)

The `verify()` card system tracks every analytical claim made in
`src/*_registry.jl`. The harness checks:

- **Inventory drift** (`test_inventory_drift.jl`) — `test/INVENTORY.jsonl`
  stays in sync with the AST of all `verify()` calls; any new or
  deleted card is caught immediately.
- **Registry coverage** (`test_atlas_logic.jl`) — every `@register`
  entry has at least one `verify()` card.
- **Documentation structure** (`test_doc_structure.jl`) — atlas hub
  pages exist for every registered quantity.

See [Identity Harness](identity-harness.md) for the full `verify()` protocol.

---

## Layer 5 — Convention lint (`test/lint/`)

Enforces the operator-convention header required in every top-level
model file (see [Conventions](../conventions.md)).

| File | Rule enforced |
|------|---------------|
| `test_convention_declarations.jl` | Every `src/models/quantum/<Model>/<Model>.jl` must contain a `# CONVENTION` block |

---

## Design Principles

1. **Every `src/` value has at least one Layer 3 card.** No analytical
   formula enters the source code without being verified against an
   independent numerical calculation.
2. **Layer 4 closes the loop.** No new `@register` entry can escape
   without a `verify()` card; no card can silently disappear.
3. **Multiple computation paths are a feature, not redundancy.** For the
   TFIM, BdG, full $2^N$ ED, and AD thermodynamics all exist as
   independent paths; agreement across all three provides strong evidence
   of correctness.
4. **Test sizes are kept small enough for exact comparison.** All
   verification tests use $N \leq 16$ for spin models, enabling
   comparison to machine precision.

---

## Further Reading

- [Cross-Verification Table](cross-checks.md) — all 8 universality ↔ model cross-checks
- [Entanglement Verification](entanglement.md) — central charge from $S(l)$
- [Disordered Systems](disordered.md) — IRFP and random singlet
- [Identity Harness](identity-harness.md) — the `verify()` card protocol
