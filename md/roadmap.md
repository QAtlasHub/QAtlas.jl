# Roadmap and Open Issues

## Architecture (v0.24 — current)

QAtlas is organised as **three parallel namespaces** over a **four-value
`status` axis**:

- `models/` — concrete Hamiltonians (`TFIM`, `XXZ1D`, …); `status` ∈
  `:exact` / `:bound` / `:approx`.
- `universalities/` — `Universality{C}` RG / CFT / RMT classes; `status =
  :universal` (derived by construction from the namespace, #659).
- `bounds/` — `Bound{D}` model-independent universal bounds (Tsirelson,
  Bekenstein, …); `status = :bound`, pinned by a `direction` and a `scheme=`
  selector when several bounds share a quantity.

A **knowledge-graph layer** (`src/core/links.jl`) makes every relation
queryable in both directions, and a **structural coherence suite**
(`src/core/coherence.jl`, C1–C9) self-reports holes: dangling references
(C1, all edge stores), kind ↔ namespace desync (C2), canonical/scheme
coherence (C3), undeclared delegation targets (C4), coverage (C6),
mutually-exclusive critical loci (C8), orphan classes (C9). Four declarative
stores feed it — `@register` (REGISTRY: how to compute a quantity),
`@realizes` (model → universality class), `@reduces` (model → model
reduction), `@about` (model cards) — all sharing one `_forward_kw_macro`
helper (#660).

**Coherence is structural, not physical.** `coherence_report()` holds
`:error` findings to **zero** (the hard invariant); `:gap` findings are
self-reported coverage holes and *need not be empty* (e.g. the `:Potts3` /
`:Potts4` orphan classes — valid 1+1D CFTs with no realizing model in the
atlas yet, #661). C5 (realization-agreement) and C7 (bound-satisfaction) are
*designed but not auto-run* slots; the physical values they would assert are
checked in `test/verification/` and `test/bounds/` via `verify` /
`verify_bound`. The coherence layer is the structural complement to those
tests, not a replacement.

### Decision: single-`Claim` unification — DEFERRED (trigger documented)

The four declarative stores could in principle collapse into one typed
`Claim` edge, with registry / fetch / verify / atlas / bibliography as views.
The v0.24 design review (#662) concluded this is **not worth doing now**, and
the evidence *strengthens* the deferral:

- The stores carry genuinely distinct invariants (REGISTRY's cross-field
  validation; REALIZES's behavioural `at` predicate + `realized_class`
  runtime mutual-exclusion). A tagged union reproduces today's per-store
  validation as per-variant validation with no net win, and weakens the
  package-load-time guards.
- The only real duplication — the four macro bodies — was already removed
  (#660) without touching the type structure.
- Adding a store is cheap: REDUCES landed without a docs parse-walker or
  `generate.jl` wiring (the docs harness is 3-way, not 4-way).

**Revisit when EITHER** (a) a fifth store lands that *also* needs cross-store
coherence checks beyond the C-suite (not expressible as a pure view like
`links.jl`), **OR** (b) REDUCES / ABOUT grow their own cross-field validation
*and* inverse-query surface. The companion `Implementation.delegates_to`
field (upgrading C4 from a presence check to an exact per-row target match,
partly subsuming C5) shares the same trigger: the first model → class
delegation on a model already carrying a `@reduces` edge.

## Remaining Issues

### Physics accuracy of documentation

The `docs/src/calc/` derivation notes provide overview-level arguments
for each result. Several notes — especially around E8 and integrable
field theory — describe **what** is true and **why** it matters, but
do not contain the full calculation. Specific gaps:

- `calc/e8-mass-spectrum-derivation.md`: the S-matrix bootstrap
  equations are stated but not solved. The actual derivation requires
  the CDD factor construction and minimal form factor program
  (Zamolodchikov 1989, Fateev 1994).
- `calc/ising-cft-magnetic-perturbation.md`: Zamolodchikov's proof
  that 8 conserved charges survive the σ perturbation is cited but not
  reproduced. This requires counting integrals of motion order by order
  in perturbation theory.
- `calc/yang-magnetization-toeplitz.md`: the Toeplitz determinant
  evaluation via Szegő's theorem is described at a high level but the
  asymptotic analysis is not shown.

These are candidates for deeper `calc/` sub-notes when expertise is
available.

### ~~Peschel method for σ^z σ^z TFIM~~  ✅ Clarified (v0.13.3)

The BdG eigenvectors for the TFIM in the σ^z σ^z convention give
**dual-fermion** correlators (Kramers-Wannier dual), not spin-basis
correlators.  For a contiguous block `A = {1, …, ℓ}` the JW unitary
factorises across the `A / Aᶜ` cut, so the spin-basis entanglement
entropy and the fermion-basis entanglement entropy are equal
(Fagotti–Calabrese 2010).  The derivation is written out in
`docs/src/calc/tfim-entanglement-peschel.md` (v0.13.3 PR #69) and the
`fetch(TFIM, VonNeumannEntropy, OBC(N); ℓ)` method uses the
correlation-matrix algorithm in O(N³) for contiguous blocks.

Full ED (`O(2^N)`) is still used by `test_entanglement_central_charge.jl`
to exercise both entropy-dispatch paths at N ≤ 16.  The duality caveat
now applies only to single-site `σ^z` correlators reached via
Pfaffians, not to contiguous-block entanglement entropy.

### ~~Central charge extraction precision~~  ✅ Tightened (v0.13.5)

PR #75 switched the tight-signal central-charge extraction to PBC
ground states via `build_tfim_sparse` / `build_spinhalf_heisenberg_sparse`
+ KrylovKit Lanczos:

| model (PBC)        | extracted c      | exact c | assertion          |
|--------------------|------------------|---------|--------------------|
| TFIM at h = J      | c(N=14) = 0.5021 | 0.5     | rtol = 0.01        |
| Heisenberg (1/N ext.) | c(∞) = 1.014 | 1.0     | rtol = 0.03        |

The OBC testset (rtol = 0.10 / 0.20) is kept for open-boundary code-
path coverage; both routes pass.

### ~~Bond counting convention for small PBC lattices~~  ✅ Decoupled (v0.13.4)

The IsingSquare transfer-matrix test no longer depends on whichever
convention `Lattice2D.bonds(lat)` adopts; PR #71 moved the bond-list
construction to `test/util/classical_partition.jl`'s
`square_pbc_bond_pairs(Lx, Ly)` and dropped the `[sources]` pin on
Lattice2D / LatticeCore.  Current `Lattice2D` (≥ 0.2.8) is the
registered release.

### ~~Two API styles~~  ✅ Resolved (v0.13)

All models have been migrated to concrete parameter-carrying structs:

- `TFIM(; J, h)`, `E8()`, `XXZ1D(; J, Δ)`, `Heisenberg1D(; J)`,
  `IsingSquare(; J, Lx, Ly)`, `Honeycomb(; t, Lx, Ly)`,
  `Kagome(; t, Lx, Ly)`, `Lieb(; t, Lx, Ly)`, `Triangular(; t, Lx, Ly)`,
  `Universality{C}`.
- Quantities are concrete structs: `Energy`, `FreeEnergy`,
  `MagnetizationX/Y/Z`, `SusceptibilityXX/YY/ZZ`, `VonNeumannEntropy`,
  `ThermalEntropy`, `RenyiEntropy(α)`, `ZZCorrelation{M}` (parametric),
  `ZZStructureFactor`, `CentralCharge`, `LuttingerParameter`,
  `LuttingerVelocity`, `FermiVelocity`, `E8Spectrum`,
  `TightBindingSpectrum`, `ExactSpectrum`, `GroundStateEnergyDensity`,
  `PartitionFunction`, `CriticalTemperature`, `CriticalExponents`,
  `GrowthExponents`, …
- `const SpinWaveVelocity = LuttingerVelocity` — type-level alias for
  the 1D critical spin chain community.
- BCs `OBC(N)`, `PBC(N)` carry their own size; `Infinite()` unchanged.
- Legacy Symbol-dispatch calls (`fetch(:TFIM, :energy, OBC(); N, J, h)`
  etc.) are routed through the deprecation shims in `src/deprecate/`
  with a `maxlog=1` info log (keyed per `(model, quantity)` pair).
  Incidental verification call sites have been ported to the concrete-
  struct API (PRs #73, #78); only the dedicated
  `"… legacy Symbol dispatch …"` testsets still reach the shims, and
  those wrap each call with `@test_logs (:info, r"symbol-dispatch") …`
  so CI output stays clean.  The directory can be `git rm`-ed
  wholesale at v1.0 cut.
- Aqua.jl static checks are enabled in `test/test_aqua.jl`
  (ambiguities, deps_compat, stale_deps, piracies).

Version bumped to 0.13.0.  See PRs #40 (foundation), #41 (TFIM + E8),
#42 (classical + tight-binding + Graphene→Honeycomb rename), #44
(XXZ1D).

## Future Directions

### Near-term (infrastructure exists)

**~~Sparse ED for larger N~~** ✅ Landed (v0.13.3+):
- `build_tfim_sparse`, `build_spinhalf_heisenberg_sparse` and
  `build_xxz_sparse` in `test/util/sparse_ed.jl`, with
  `ground_state_krylov` driving `KrylovKit.eigsolve`.  The PBC
  entanglement-entropy, TFIM-BdG-vs-ED, z-exponent, and Luttinger-K
  cross-checks now run at N ∈ {14, 16} on the default CI profile.
  N = 18–20 remain accessible under `QATLAS_TEST_FULL = 1` (nightly).

**More universality classes from Wikipedia table**:
- Directed percolation (d=1,2,3: numerical, d≥4: MF)
- Self-avoiding walk (d=2: exact, d=3: numerical)
- Conserved directed percolation (Manna class)
- Tricritical Ising (c = 7/10)
- These are all data-entry tasks with the `Universality{C}` framework.

**Finite-size scaling util** (`test/util/fss.jl`):
- Systematic extraction of critical exponents from ED data at multiple
  N values. Binder cumulant crossing, data collapse, etc.
- Would make the universality cross-checks more rigorous.

### Medium-term (requires new packages or significant work)

**ITensors.jl DMRG benchmarks**:
- Compare DMRG ground state energy against QAtlas exact values
  (Heisenberg Bethe ansatz, TFIM BdG).
- Requires ITensors.jl as a test dependency.
- Goal: QAtlas as the "ground truth" for tensor network benchmarks.

**Disordered systems — quantitative IRFP**:
- Extract c_eff = ln 2 from the random TFIM at the Fisher IRFP
  with proper disorder averaging (100+ realizations, N ~ 14).
- Random singlet phase of the Heisenberg chain: verify the
  log-averaged correlation function exponent.
- Harris criterion test: compare clean α to disorder relevance.

**Anderson localization** (1D TB + random on-site potential):
- All states localized in 1D (exact result).
- Localization length ξ ∝ 1/W² for weak disorder.
- Can use existing `build_tight_binding` + random diagonal.

**Entanglement entropy for free fermions** (Peschel done right):
- Implement the correlation-matrix method for the σ^x σ^x convention
  TFIM (direct JW, no KW duality needed).
- O(N³) computation enables N ~ 1000 for precise c extraction.

### Long-term (research-level)

**TFIM + longitudinal field (TFIML)**:
- Numerical verification of E8 mass ratios on finite chains.
- Requires computing the excitation spectrum (not just ground state)
  of the TFIML Hamiltonian and identifying the 8 stable particles.
- Connection to `src/universalities/E8.jl`.

**2D models on Lattice2D**:
- Classical Monte Carlo verification of Onsager T_c, Yang M(T) via
  Lattice2DMonteCarlo.
- Finite-size scaling of the Binder cumulant crossing → extract ν.
- Requires fixing the Lattice2DMonteCarlo LatticeCore compat issue.

**Bethe ansatz beyond ground state**:
- XXZ chain Luttinger parameters (K, u) — **done in v0.13** for the
  full critical regime −1 < Δ ≤ 1.  Ground-state energy density
  currently exposes only the three exact points Δ ∈ {−1, 0, 1}.
- General-Δ Yang-Yang ground-state integral (Δ ∈ (−1, 1) smoothly
  between the exact points) — deferred to a follow-up PR; several
  equivalent forms in the literature differ by nontrivial
  substitutions and sign conventions, requires careful boundary-limit
  verification.
- Heisenberg: des Cloizeaux-Pearson dispersion ε(k) = (π/2)|sin k| —
  needed for dynamical structure factor; `LuttingerVelocity` already
  returns the right Δ = 1 limit.
- Finite-temperature Bethe ansatz (TBA) for thermodynamic quantities.

**Kitaev models**:
- Kitaev chain (1D): topological phase boundary, Majorana edge modes.
  Same BdG framework as TFIM.
- Kitaev honeycomb (2D): exact Majorana solution, A/B phase diagram.
  Requires 2D lattice Majorana implementation.

**Detailed E8 derivation**:
- Fill in the missing calculation steps in `calc/e8-mass-spectrum-derivation.md`:
  bootstrap equation solution, CDD factors, minimal form factors.
- Requires expertise in integrable field theory.

## Development Practices

- **Version bump**: patch per PR, minor for significant features.
- **Testing**: 3-tier (standalone, verification, cross-verification).
  Every `src/` value must have independent test verification.
- **Documentation**: Zettelkasten `calc/` notes for derivations, linked
  from model/method/universality pages.
- **Citations**: each stored value traces to a specific paper + equation.
- **BLAS threading**: `runtests.jl` sets `BLAS.set_num_threads(Sys.CPU_THREADS)`.
- **Julia threads**: invoke `Pkg.test` with
  `julia_args=\`-t auto --heap-size-hint=96G\`` on the 36-core / 128GB
  host to avoid swap thrash; default 1-thread runs OOM on the
  entanglement-central-charge tests.
- **Adaptive sizing**: ED tests use larger N when RAM > 50GB.
