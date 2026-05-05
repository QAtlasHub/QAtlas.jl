# Universality Classes

## Overview

At a continuous phase transition, thermodynamic quantities diverge
(or vanish) as power laws characterised by **critical exponents**.
A universality class is the set of all systems that share the same
critical exponents. Membership depends only on:

1. **Spatial dimension** $d$
2. **Symmetry** of the order parameter
3. **Range of interactions** (short-range vs long-range)

Microscopic details --- lattice structure, coupling strength, even
the distinction between classical and quantum --- are irrelevant.
This remarkable fact, grounded in the renormalization group, is the
reason a single library entry can serve as a reference for an entire
family of physical models.

---

## Universality Classes in QAtlas

| Class | Symbol | Dimensions | Type | Key feature | Page |
|-------|--------|-----------|------|-------------|------|
| Ising | `:Ising` | $d = 2, 3, \geq 4$ | Exact / Bootstrap / MF | $\mathbb{Z}_2$ symmetry | [→](ising.md) |
| Percolation | `:Percolation` | $d = 2, 3, \geq 6$ | Exact / MC / MF | Geometric transition | [→](percolation.md) |
| 3-state Potts | `:Potts3` | $d = 2$ | Exact | $S_3$ symmetry | [→](potts.md) |
| 4-state Potts | `:Potts4` | $d = 2$ | Exact | Marginal, log corrections | [→](potts.md) |
| KPZ | `:KPZ` | $1+1$D | Exact | Non-equilibrium growth | [→](kpz.md) |
| XY | `:XY` | $d = 2, 3, \geq 4$ | BKT / Bootstrap / MF | $\mathrm{O}(2)$ symmetry | [→](on-models.md) |
| Heisenberg | `:Heisenberg` | $d = 3, \geq 4$ | Bootstrap / MF | $\mathrm{O}(3)$ symmetry | [→](on-models.md) |
| Mean-Field | `MeanField()` | $d \geq d_c$ | Exact | Baseline reference | [→](mean-field.md) |
| E8 | `:E8` | — | Exact mass ratios | Integrable field theory | [→](e8.md) |

---

## The `Universality{C}` API

All universality data in QAtlas is accessed through the type-safe
`Universality{C}` interface:

```julia
using QAtlas

# Construct a universality-class handle
u = Universality(:Ising)   # equivalent to Universality{:Ising}()

# Fetch critical exponents for a given dimension
e = QAtlas.fetch(u, CriticalExponents(); d=2)
# (β = 1//8, ν = 1//1, γ = 7//4, η = 1//4, δ = 15//1, α = 0//1, c = 1//2)

# Mean-field has its own type (no dimension needed)
e_mf = QAtlas.fetch(MeanField(), CriticalExponents())
```

Return types depend on whether the exponents are analytically known:

| Situation | Return type | Example |
|-----------|-------------|---------|
| Exact ($d = 2$ Ising, Potts, ...) | `Rational{Int}` | `β = 1//8` |
| Numerical ($d = 3$ bootstrap, MC) | `Float64` + `_err` fields | `β = 0.32642, β_err = 1e-5` |
| Mean-field ($d \geq d_c$) | `Rational{Int}` | `β = 1//2` |

---

## Scaling Relations

For any universality class with standard critical exponents, the
following relations hold (see [calculation note](../calc/ising-scaling-relations.md)):

$$\alpha + 2\beta + \gamma = 2 \qquad \text{(Rushbrooke)}$$
$$\gamma = \beta(\delta - 1) \qquad \text{(Widom)}$$
$$\gamma = \nu(2 - \eta) \qquad \text{(Fisher)}$$
$$2 - \alpha = d\nu \qquad \text{(Josephson hyperscaling, } d < d_c\text{)}$$

QAtlas stores exact-class exponents as `Rational{Int}`, so these
identities can be verified with **zero floating-point error** in
the test suite.

---

## Further Reading

- [Ising](ising.md) --- the most thoroughly studied class
- [Mean-Field](mean-field.md) --- baseline for $d \geq d_c$
- [Cross-Verification](../verification/cross-checks.md) --- how
  universality predictions are tested against model-specific results
