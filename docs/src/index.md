# QAtlas.jl

**QAtlas** (QUAntum Reference Table for Exact Tests) is a curated
dictionary of rigorous results in quantum and statistical physics.
Every stored value is traced to a specific publication and
cross-validated against independent calculations.

## What Makes QAtlas Different

Unlike typical numerical libraries, QAtlas focuses on **authoritative
reference values** — exact analytical results, high-precision
conformal bootstrap bounds, and Bethe ansatz solutions. Each value
is accompanied by:

1. **Precise citation**: author, year, journal, equation number
2. **Derivation sketch**: enough to independently verify
3. **Cross-validation**: tested against independent computation
4. **Connections**: linked to universality classes and other models

## Quick Start

```julia
using QAtlas

# Onsager critical temperature
Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature())

# TFIM ground-state energy
E₀ = QAtlas.fetch(:TFIM, :energy, OBC(); N=16, J=1.0, h=0.5)

# 2D Ising universality: exact exponents (Rational)
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
# (β = 1//8, ν = 1//1, γ = 7//4, η = 1//4, ...)
```

## Contents

### Models

Exact solutions for specific physical models.

| Model | Type | Key Results | Page |
|-------|------|-------------|------|
| TFIM | Quantum | Energy, gap, thermal observables, entanglement | [→](models/quantum/tfim.md) |
| IsingSquare | Classical | $Z$, $T_c$, $M(T)$ | [→](models/classical/ising-square.md) |
| Heisenberg1D | Quantum | Dimer, 4-site PBC, Bethe $e_0$ | [→](models/quantum/heisenberg.md) |
| MajumdarGhosh | Quantum | Exact dimer GS at $J_2/J_1=1/2$, $E_0/N=-3J/8$ | [→](models/quantum/majumdar_ghosh.md) |
| XXZ1D | Quantum | Exact $\Delta \in \{-1, 0, 1\}$ + Luttinger $K, u$ | [→](models/quantum/xxz.md) |
| Honeycomb TB | Quantum | Bloch spectrum (honeycomb / graphene) | [→](models/quantum/tightbinding/honeycomb.md) |
| Kagome TB | Quantum | Flat band at $+2t$ | [→](models/quantum/tightbinding/kagome.md) |
| Lieb TB | Quantum | Flat band at $E = 0$ | [→](models/quantum/tightbinding/lieb.md) |
| Triangular TB | Quantum | Frustrated band $[-6t, +3t]$ | [→](models/quantum/tightbinding/triangular.md) |

### Universality Classes

Critical exponents and scaling relations via `Universality{C}`.

| Class | Dimensions | Type | Page |
|-------|-----------|------|------|
| Ising | $d = 2, 3, \geq 4$ | Exact / Bootstrap / MF | [→](universalities/ising.md) |
| Percolation | $d = 2, 3, \geq 6$ | Exact / MC / MF | [→](universalities/percolation.md) |
| Potts ($q = 3, 4$) | $d = 2$ | Exact | [→](universalities/potts.md) |
| KPZ | $1+1$D | Exact | [→](universalities/kpz.md) |
| XY / Heisenberg | $d = 2, 3, \geq 4$ | BKT / Bootstrap / MF | [→](universalities/on-models.md) |
| E8 | — | Exact mass ratios | [→](universalities/e8.md) |

### Verification

How QAtlas ensures physical correctness.

- [Philosophy](verification/index.md) — three-layer testing strategy
- [Cross-Checks](verification/cross-checks.md) — universality ↔ model connections
- [Entanglement](verification/entanglement.md) — central charge from $S(l)$
- [Disordered Systems](verification/disordered.md) — IRFP, random singlet

### Methods

Computational techniques used by QAtlas, with physical justification.

- [Transfer Matrix](methods/transfer-matrix/index.md)
- [Bloch Hamiltonian](methods/bloch-hamiltonian/index.md)
- [Exact Diagonalization](methods/exact-diagonalization/index.md)
- [Automatic Differentiation](methods/automatic-differentiation/index.md)
- [Calabrese-Cardy Formula](methods/calabrese-cardy/index.md)

### API Reference

```@autodocs
Modules = [QAtlas]
```
