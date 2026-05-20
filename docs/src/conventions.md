# Operator Conventions

This page is the **normative specification** for what numerical values
`fetch(model, ::AbstractQuantity, bc; …)` is allowed to return in
QAtlas.jl.  Every model file under `src/models/quantum/` must declare a
`CONVENTION` header that conforms to this policy (CI enforces this — see
`test/lint/test_convention_declarations.jl`).

The motivation is reader unambiguity and downstream interoperability:
the same `Quantity` type returned by two different hubs must mean the
same physical thing, in the same units, with the same operator
normalisation.  When that promise breaks, `verify()` cards comparing the
two hubs silently fail (or worse, silently pass at the wrong value).

## Policy summary

| Model family | Observable convention |
|---|---|
| Spin / qubit (Heisenberg, XXZ, TFIM, Cluster, Kitaev, Compass, Toric, …) | **Spin operator** `S^α` with eigenvalues in `[-S, +S]` |
| Fermionic (Hubbard, Schwinger, tight-binding, Gross-Neveu, SYK, …) | **Number / bilinear** `n_i = c†_i c_i`, `⟨c†_i c_j⟩` |
| Hard-core boson / Rydberg / constrained (PXP, …) | **Occupation** `n_i ∈ {0, 1}` for density family; spin S for `MagnetizationX/Y/Z`-family |
| Topological / anyonic (Toric, XCube, ChernSimons, FibonacciAnyons, …) | **Operator-product** expectation (Wilson loop, GSD, S-matrix entries, …) |

These branches do not overlap in their natural observable categories,
so there is no ambiguity in practice.  When a model spans two families
(e.g. PXP is both "spin Hamiltonian written with σ" and "Rydberg
density observables"), the **observable** category wins — the header
declares the convention for `fetch` return values, not the Hamiltonian
formula.

## Spin convention

For spin systems QAtlas standardises on the **spin operator** convention.
For spin-S systems the spin operators `S^x, S^y, S^z` satisfy
`[S^α, S^β] = iε^{αβγ} S^γ` and have eigenvalues in `{-S, -S+1, …, +S}`.
For spin-1/2 systems, `S^α = σ^α / 2`.

This is the convention you find in modern condensed-matter textbooks
(Sachdev, Auerbach, Giamarchi) and it generalises uniformly to S=1, S=3/2, …

### Observable scaling table (σ → S, spin-1/2)

`n` is the number of spin operators in the observable.  The S-convention
value equals the σ-convention value times `2^{-n}`.

| Quantity (`AbstractQuantity` subtype) | `n` | Conversion factor | Example (spin-1/2, T → ∞ saturation) |
|---|---|---|---|
| `MagnetizationX`, `MagnetizationY`, `MagnetizationZ`, `…Local` | 1 | `× 1/2` | σ-conv `1` → S-conv `1/2` |
| `XXCorrelation`, `YYCorrelation`, `ZZCorrelation` | 2 | `× 1/4` | σ-conv `1` (i = j) → S-conv `1/4` |
| `XXStructureFactor`, `YYStructureFactor`, `ZZStructureFactor` | 2 | `× 1/4` | σ-conv `1` → S-conv `1/4` |
| `SusceptibilityXX`, `SusceptibilityYY`, `SusceptibilityZZ` | 2 | `× 1/4` | σ-conv `β/N` (Curie) → S-conv `β/(4N)` |
| `StringOrderParameter` (when defined via spin-1/2 σ) | `k` | `× 2^{-k}` | depends on string length |

For higher-spin systems (S=1, S=3/2, …) the σ representation does not
exist; `fetch` is unambiguous because there is only one option.  The
diagonal correlator at infinite temperature, by SU(2) symmetry, is
`⟨S^α_i S^α_i⟩_∞ = S(S+1)/3` per site (e.g. `2/3` for S=1; `1/4` for S=1/2).

### Quantities NOT affected by σ vs S choice

The following return values are independent of any spin operator
normalisation and so the convention discussion does not apply to them.

`Energy`, `FreeEnergy`, `SpecificHeat`, `MassGap`, `ChargeGap`,
`SpinGap`, `CorrelationLength`, `ResidualEntropy`,
`ThermalEntropy`, `VonNeumannEntropy`, `RenyiEntropy`,
`FidelitySusceptibility`, `LoschmidtEcho`, `CentralCharge`,
`ConformalWeights`, `PrimaryFields`, `E8Spectrum`, `FractalDimension`,
`LuttingerParameter`, `FermiVelocity`, `LuttingerVelocity`,
`TopologicalInvariant`, `EdgeModeEnergy`, `GroundStateDegeneracy`,
`TopologicalEntanglementEntropy`, `AnyonStatistics`,
`WignerSurmise`, `TracyWidom`, `MeanRatio`, `SpectralFormFactor`,
`CasimirEnergyCorrection`, `SteadyStateCurrent`, `ChiralCondensate`.

The `Energy*` family does carry a *coupling-constant* convention (the
Hamiltonian in σ vs S has different J, h numerics) but the return value
of `fetch(model, Energy(...), bc)` is just a number with units of `[J]`
or `[h]` — the model's struct field values define the units, not a
universal QAtlas-wide rule.

### Hamiltonian vs observable

The convention **for `fetch` return values** is fixed per QAtlas-wide
policy.  The convention **inside the Hamiltonian docstring** is the
choice of the model author — typically the one literature uses for that
model.  TFIM keeps its Pfeuty σ Hamiltonian; Heisenberg / XXZ keep
their `S · S` form.  The header block must declare both:

```julia
# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide policy)
```

This is so a reader who lands on `TFIM.jl` and sees Pfeuty's
`m_x = (1 − (h/J)^2)^{1/8}` in a docstring does not silently apply that
formula as the value returned by `fetch(TFIM(...), MagnetizationX(),
…)` — which under this policy is half that.

## Fermion convention

For models whose elementary degrees of freedom are fermions (`Hubbard1D`,
`ExtendedHubbard1D`, `TightBinding1D`, `TightBindingV1D`,
`SchwingerModel`, `GrossNeveu`, `SYK`, …), observables are defined in
terms of fermion bilinears:

| Quantity | Definition |
|---|---|
| Number / density | `n_i = c†_i c_i` (or `n_{i,σ} = c†_{i,σ} c_{i,σ}` per spin flavour) |
| Spin density (if spin index exists) | `S^z_i = (n_{i,↑} − n_{i,↓})/2`, `S^+_i = c†_{i,↑} c_{i,↓}` |
| Pair amplitude | `Δ_i = c_{i,↓} c_{i,↑}` |
| Hopping correlator | `⟨c†_i c_j⟩` |
| Chiral condensate | `⟨ψ̄ψ⟩` (Dirac convention) |

Note: derived spin observables `MagnetizationX/Y/Z` on fermionic models
follow the **spin S convention** (factor `1/2` already baked into the
definition), so the spin convention table above applies uniformly.

## Hard-core boson / Rydberg / constrained spin (PXP)

Some models have a spin-1/2 Hamiltonian formally written in σ-operators
but whose natural observables are *occupations* `n_i ∈ {0, 1}` rather
than spin projections.  The canonical example is **PXP1D** (Rydberg
blockade chain):

```
H = Ω Σ_i P_{i-1} σ^x_i P_{i+1},  P_i = (1 − σ^z_i)/2,  n_i = (1 + σ^z_i)/2
```

For these models, the Rydberg literature reports `⟨n_i⟩` (excited-atom
density), `⟨n_i n_j⟩` (density-density correlator), Z₂ order parameter,
etc.  Forcing them into "S^z = ⟨n⟩ − 1/2" hurts readability against the
literature without buying any unification gain (PXP-specific quantities
will never be cross-checked against a Heisenberg solver).

**Policy for these models:**

1. If a quantity is a spin observable in the formal sense
   (`MagnetizationX/Y/Z`, `Susceptibility…`, `*Correlation`,
   `*StructureFactor`) it follows the **spin S convention** as for any
   other spin system.  This is so cross-model `verify()` cards can still
   compare against an S=1/2 reference.
2. For literature-canonical Rydberg observables (Rydberg density,
   density-density correlator, Néel order parameter, scar revival
   probability), the model file may introduce occupation-convention
   quantities and document them as `n ∈ [0, 1]`.  These are model-family
   quantities, not the universal-spin ones.

In practice, PXP1D's Phase 1 (`Energy{:per_site}` only) and Phase 2 (scar
diagnostics) sit cleanly outside the spin-observable surface, so the
choice does not bite the current API.  The Phase 2 scar Z₂ survival
probability is a *projector* expectation `⟨ψ(t) | Z_2 ⟩⟨ Z_2 | ψ(t)⟩`,
which is a probability `∈ [0,1]` independent of any spin/occupation
convention question.

## Topological / operator-product

For topological models the observables of interest are typically *not*
local spin expectations:

- `GroundStateDegeneracy` is a positive integer (the dimension of the
  topological ground-state subspace) — no convention needed.
- `TopologicalEntanglementEntropy` is a real number (a difference of
  von Neumann entropies) — no convention needed.
- `AnyonStatistics` is a complex phase or matrix (R/F symbols) — defined
  by the anyon model itself.
- Wilson loop / string operator expectations are products of stabilizer
  generators — the value is the eigenvalue of the operator product on
  the state, already convention-free.

Stabilizer operators (Toric Code's `A_v = ∏ σ^x`, `B_p = ∏ σ^z`;
ChernSimons modular data; Kitaev plaquette operators) are themselves
*products* of σ-matrices with eigenvalues `±1`.  Their expectation
values are dimensionless and do not pick up any factor under σ ↔ S
rescaling because each σ inside the stabilizer is paired with a σ in
the *same* operator (not a separate observable).

## Required header in every model file

Every file `src/models/quantum/<Model>/<Model>.jl` (i.e. the top-level
model definition file, not the per-quantity solver files) must contain
a `CONVENTION` block exactly matching this regex:

```
^#\s*CONVENTION\s*$
(^#.*$\n){2,}
```

i.e. a comment line saying `# CONVENTION`, followed by at least two more
comment lines naming the Hamiltonian convention and the observable
convention.  Example accepted forms:

```julia
# CONVENTION
#   Hamiltonian: Pauli σ (this file)
#   Observable:  Spin S = σ/2  (QAtlas-wide spin convention)
```

```julia
# CONVENTION
#   Hamiltonian: Spin S (this file)
#   Observable:  Spin S        (QAtlas-wide spin convention)
```

```julia
# CONVENTION
#   Hamiltonian: Pauli σ with Rydberg projectors (this file)
#   Observable:  Spin S for MagnetizationX/Y/Z-family,  occupation n ∈ [0,1] for density-family
#   Reference:   docs/src/conventions.md §Hard-core boson / Rydberg
```

```julia
# CONVENTION
#   Hamiltonian: Fermion bilinears c†c
#   Observable:  Fermion conventions (number n = c†c, bilinear ⟨c†c⟩); derived
#                spin observables MagnetizationX/Y/Z follow spin S = σ/2
```

```julia
# CONVENTION
#   Hamiltonian: Stabilizer (A_v, B_p products of Pauli σ)
#   Observable:  Operator-product expectations (Wilson loops, GSD, TEE);
#                no spin/occupation choice required
```

If a file does not need its own convention block because it only
contains per-quantity solver methods (e.g. `TFIM_dynamics.jl`,
`XXZ_thermal.jl`), it is exempt — only the top-level model file (the
one defining the `struct <Model> <: AbstractQAtlasModel`) must declare
the convention.

## Migration of existing σ-convention observables

The following hubs currently return σ-convention values for one or more
observables, and **must** be migrated to S convention.  Each migration
is a separate PR with:

- Solver implementation updated to apply the `2^{-n}` factor
- Every affected `verify()` card's `independent` value updated
- Refs strings updated to note the S-convention return value
- Header block added or updated

| Hub | Affected observables | PR slot |
|---|---|---|
| `TFIM` | `MagnetizationX(Local)`, `*StructureFactor`, `*Correlation`, `Susceptibility*` | open |
| `Cluster1D` | (mainly Energy / MassGap; observables limited) | tbd |
| `Compass1D` | (mainly MassGap; observables limited) | tbd |
| `LongRangeIsing1D`, `MixedFieldIsing1D`, `XYh1D`, `LongRangeXY1D` | spin observables | tbd |
| `Kitaev1D`, `KitaevHoneycomb` | spin observables (if any beyond Energy / MassGap) | tbd |
| `ToricCode`, `XCube` | (mainly GSD / TEE; convention-free) | none required |

After the migration cascade, every `verify()` card whose `independent`
value would change must be re-pushed.  CI will catch any divergence
because the existing cards exercise both the analytical formula and
the solver value.
