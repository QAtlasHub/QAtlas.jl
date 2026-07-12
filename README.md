# QAtlas.jl

> ⚠️ **AI-assisted draft — feedback welcome.**
> The source code and the derivation notes under `docs/src/` are written
> with heavy LLM assistance (Claude) and are being independently
> reviewed. Found an error — a formula, a proof step, a citation, a
> type signature? Please
> [**open an issue**](https://github.com/QAtlasHub/QAtlas.jl/issues/new?labels=docs&title=%5Bdocs%5D%20error%20report)
> or a PR. Recent review cycles have already caught one flagship bug
> (E₈ mass ratios `m₇`, `m₈` were 2× the literature value, fixed in
> [PR #83](https://github.com/QAtlasHub/QAtlas.jl/pull/83));
> more eyes keep making the package better.

[![docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://codes.sota-shimozono.com/QAtlas.jl/stable/)
[![docs: dev](https://img.shields.io/badge/docs-dev-purple.svg)](https://codes.sota-shimozono.com/QAtlas.jl/dev/)
[![Julia](https://img.shields.io/badge/julia-v1.12+-9558b2.svg)](https://julialang.org)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[![codecov](https://codecov.io/gh/QAtlasHub/QAtlas.jl/graph/badge.svg?token=dnvvwd7DVo)](https://codecov.io/gh/QAtlasHub/QAtlas.jl)
[![Build Status](https://github.com/QAtlasHub/QAtlas.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/QAtlasHub/QAtlas.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**QAtlas** (QUAntum Reference Table for Exact Tests) is a Julia
package exposing a curated dictionary of **rigorous results** in
quantum and statistical physics. Every stored value traces back to a
specific publication and is cross-validated against independent
numerical computations (dense / sparse exact diagonalisation,
brute-force enumeration, Bloch diagonalisation, automatic
differentiation, bipartite fluctuations, …).

> **Confidence tiers.**
> Results the maintainer can independently derive and has checked
> line-by-line (TFIM BdG + thermodynamics, Onsager $T_c$, Yang $M(T)$,
> Heisenberg dimer, bipartite-fluctuation Luttinger $K$, PBC
> Calabrese–Cardy central charge) carry the highest confidence.
> Others — tight-binding Bloch formulas, 3D O(N) bootstrap values,
> $E_8$ mass ratios — are cross-checked against independent numerical
> computations but the underlying derivations have not all been
> re-done by hand.
> **If you use QAtlas values in a publication, please verify them
> against the cited original references.**

## Quick Start

```julia
using QAtlas

# ── Classical Ising on the 2D square lattice ─────────────────────
# Onsager (1944) — critical temperature
Tc = QAtlas.fetch(IsingSquare(; J = 1.0), CriticalTemperature())
#  = 2J / ln(1 + √2)

# Yang (1952) — spontaneous magnetisation
M = QAtlas.fetch(IsingSquare(; J = 1.0), SpontaneousMagnetization(); β = 0.5)

# Transfer-matrix partition function on an Lx × Ly torus
Z = QAtlas.fetch(IsingSquare(; J = 1.0, Lx = 4, Ly = 4), PartitionFunction(); β = 0.44)

# ── Universality classes ─────────────────────────────────────────
# 2D Ising — exact rational exponents
e2 = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d = 2)
#  (α = 0//1, β = 1//8, γ = 7//4, δ = 15//1, ν = 1//1, η = 1//4, c = 1//2)

# 3D Ising — Kos-Poland-Simmons-Duffin-Vichi 2016 bootstrap
e3 = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d = 3)
#  (α = 0.11009, α_err = 1e-5, β = 0.32642, β_err = 1e-5, …)

# E₈ mass spectrum (Zamolodchikov 1989) — normalised m₁ = 1
masses = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
#  [1.0, φ, 1.989, 2.405, 2.956, 3.218, 3.891, 4.783]

# ── Quantum spin chains ──────────────────────────────────────────
# Hulthén (1938) — AF Heisenberg ground-state energy density
e0 = QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity())  # J·(1/4 − ln 2)

# TFIM BdG thermal observables (Pfeuty 1970 + free-fermion thermo)
ε = QAtlas.fetch(TFIM(; J = 1.0, h = 0.5), Energy(), OBC(24); beta = 5.0)
Δ = QAtlas.fetch(TFIM(; J = 1.0, h = 2.0), MassGap(), Infinite())   # = 2|h−J|

# Bipartite-block von Neumann entropy (Peschel method, OBC)
S = QAtlas.fetch(TFIM(; J = 1.0, h = 1.0), VonNeumannEntropy(), OBC(16); ℓ = 8)

# XXZ chain — Luttinger parameter + spin-wave velocity (|Δ| < 1)
K = QAtlas.fetch(XXZ1D(; J = 1.0, Δ = 0.5), LuttingerParameter(), Infinite())  # = 3/4
u = QAtlas.fetch(XXZ1D(; J = 1.0, Δ = 1.0), SpinWaveVelocity(), Infinite())    # = π/2

# ── Tight-binding band structures ────────────────────────────────
# Honeycomb (Graphene) — exported alias `Graphene === QAtlas.Honeycomb`
λ = QAtlas.fetch(Graphene(; t = 1.0, Lx = 6, Ly = 6), TightBindingSpectrum())
# Kagome / Lieb / Triangular need the QAtlas. prefix (Lattice2D name clash)
λ_kagome = QAtlas.fetch(QAtlas.Kagome(; t = 1.0, Lx = 6, Ly = 6), TightBindingSpectrum())
```

## What's Inside

### Models

| Model | Example quantities | Method |
| ----- | ---------- | ------ |
| **`IsingSquare`** | `PartitionFunction`, `CriticalTemperature`, `SpontaneousMagnetization` | Transfer matrix · Onsager 1944 · Yang 1952 |
| **`TFIM`** | `Energy`, `FreeEnergy`, `ThermalEntropy`, `SpecificHeat`, `MassGap`, `MagnetizationX`, `SusceptibilityXX/ZZ`, `ZZCorrelation{:static,:dynamic,:lightcone}`, `VonNeumannEntropy` | Jordan-Wigner + BdG (Pfeuty 1970) · Peschel correlation matrix |
| **`XXZ1D`** | `Energy` (Δ ∈ {−1, 0, 1}), `LuttingerParameter`, `LuttingerVelocity` / `SpinWaveVelocity`, `FermiVelocity`, `CentralCharge` | Bethe ansatz closed forms · bosonisation |
| **`Heisenberg1D`** | `ExactSpectrum` (N=2 OBC, N=4 PBC), `GroundStateEnergyDensity` | Bethe ansatz (Hulthén 1938) |
| **`Graphene`** / **`QAtlas.Honeycomb`** | `TightBindingSpectrum` | Bloch Hamiltonian (honeycomb) |
| **`QAtlas.Kagome`** | Spectrum + flat band at +2t | Bloch 3×3 |
| **`QAtlas.Lieb`** | Spectrum + flat band at E=0 | Bloch 3×3 |
| **`QAtlas.Triangular`** | Spectrum (frustrated) | Bloch scalar |
| **`E8`** | `E8Spectrum` — 8 mass ratios | Perron-Frobenius of $E_8$ Cartan matrix (Zamolodchikov 1989) |

Additional lattices (**Dice**, **UnionJack**, **ShastrySutherland**) are covered in the verification suite via the generic Bloch builder in `src/util/bloch.jl`.

### Universality classes

Access via `Universality{C}` at dimension `d`:

| Class | d | Type | Key exponents |
| ----- | - | ---- | ------------- |
| **Ising** | 2 | Exact (Rational) | β = 1/8, ν = 1, c = 1/2 |
| | 3 | Bootstrap (Float64 + err) | β ≈ 0.32642(1), ν ≈ 0.62997(1) |
| | ≥ 4 | Mean-field | β = 1/2, ν = 1/2 |
| **Percolation** | 2 | Exact | β = 5/36, ν = 4/3 |
| | 3 | Numerical | |
| **3-state Potts** | 2 | Exact | β = 1/9, ν = 5/6 |
| **4-state Potts** | 2 | Exact | β = 1/12, ν = 2/3 |
| **XY** | 2 | BKT | η = 1/4 (universal jump) |
| | 3 | Bootstrap (O(2)) | β ≈ 0.34869(7) |
| **Heisenberg** | 3 | Bootstrap (O(3)) | β ≈ 0.3689(3) |
| **KPZ** | 1 + 1 | Exact | z = 3/2, β = 1/3, α = 1/2 |
| **Mean-field** | any | Exact | β = 1/2, ν = 1/2 |
| **E₈** | — | Exact (integrable FT) | 8 mass ratios |

Exact exponents are `Rational{Int}`, so scaling relations hold **exactly** (no floating-point error):

```julia
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d = 2)
e.α + 2 * e.β + e.γ == 2        # Rushbrooke
e.γ == e.β * (e.δ - 1)           # Widom
e.γ == e.ν * (2 - e.η)           # Fisher
2 - e.α == 2 * e.ν               # Josephson
```

## Verification

QAtlas is **self-verifying**: every stored value in `src/` is cross-checked against independent computations in `test/`. At v0.13.6 the suite runs **1124 tests in ~2 min 30 s** on a 4-core CI runner (nightly with `QATLAS_TEST_FULL=1` extends to N = 16–20 sparse ED).

Layers of verification:

- **Standalone** — special values, limits, and scaling-relation identities asserted against closed-form rationals.
- **Literature cross-check** (`test/verification/test_universality_literature_values.jl`) — numerical universality exponents match the cited paper's tabulated decimals at the quoted precision (Kos 2016, Chester 2020/2021).
- **Lattice ED verification** — real-space exact diagonalisation (dense or sparse + KrylovKit Lanczos) on [Lattice2D.jl](https://github.com/sotashimozono/Lattice2D.jl) geometries for every tight-binding and TFIM / Heisenberg / XXZ / E₈ result.
- **AD verification** — ForwardDiff-based derivation of thermodynamic quantities from `log Z`, cross-checked against direct ensemble averages.
- **Cross-universality** — β = 1/8 extracted from Yang $M(T)$, z = 1 from the TFIM gap (PBC sparse ED, `rtol = 0.01`), c = 1/2 from PBC entanglement-entropy scaling, Luttinger $K$ from bipartite-fluctuation ED — each route tests a *different* stored dictionary entry.

## ForwardDiff compatibility

Every $Z$- and $\log Z$-like quantity accepts `ForwardDiff.Dual` numbers, so thermodynamic quantities can be obtained by differentiation. `ForwardDiff` is not a runtime dependency of QAtlas — install it alongside if you want this:

```julia
using Pkg; Pkg.add("ForwardDiff")
using ForwardDiff, QAtlas

logZ(β) = log(QAtlas.fetch(IsingSquare(; J = 1.0, Lx = 3, Ly = 3),
                           PartitionFunction(); β = β))

β0 = 0.5

# Internal energy ⟨E⟩ = -∂(log Z)/∂β
E_avg = -ForwardDiff.derivative(logZ, β0)

# Heat capacity  C_v = β² · ∂²(log Z)/∂β²
d2 = ForwardDiff.derivative(β -> ForwardDiff.derivative(logZ, β), β0)
Cv = β0^2 * d2
```

The full fluctuation-dissipation derivation from $\ln Z$ to $(F, \langle E\rangle, C_v, S)$ is written out in [`docs/src/calc/ad-thermodynamics-from-z.md`](docs/src/calc/ad-thermodynamics-from-z.md).

## Installation

```julia
using Pkg
Pkg.add("QAtlas")
```

Requires Julia ≥ 1.12. Documentation: [codes.sota-shimozono.com/QAtlas.jl/stable](https://codes.sota-shimozono.com/QAtlas.jl/stable/).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the depth standard applied to derivations in `docs/src/calc/`, the citation format, and the concrete-struct API conventions. Corrections to physics content, scaling-relation cross-checks, and independent-source confirmations are especially welcome.

## License

MIT. See [LICENSE](LICENSE).
