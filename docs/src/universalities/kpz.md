# KPZ Universality Class

## Overview

The Kardar-Parisi-Zhang (KPZ) universality class describes
**non-equilibrium stochastic growth** of interfaces. It is
fundamentally different from equilibrium critical phenomena: there is
no partition function, no free energy, and the relevant exponents
characterise dynamic scaling of a growing surface rather than
static thermodynamic singularities.

The KPZ equation in $d$ spatial dimensions is

$$\frac{\partial h}{\partial t} = \nu_0 \nabla^2 h + \frac{\lambda}{2}(\nabla h)^2 + \eta(\mathbf{x}, t)$$

where $h(\mathbf{x}, t)$ is the interface height, $\nu_0$ is a
smoothing coefficient, $\lambda$ is the non-linear growth coupling,
and $\eta$ is Gaussian white noise with
$\langle\eta(\mathbf{x},t)\eta(\mathbf{x}',t')\rangle = 2D\,\delta^d(\mathbf{x}-\mathbf{x}')\delta(t-t')$.

**Systems in this class**: ballistic deposition, Eden growth,
polynuclear growth, directed polymers in random media, TASEP
(totally asymmetric simple exclusion process).

---

## $1+1$ Dimensions --- Exact Exponents

In $1+1$D (one spatial + one temporal dimension), the KPZ exponents
are known exactly.

### Growth Exponents

| Exponent | Value | Definition | Reference |
|----------|-------|-----------|-----------|
| $\beta$ | $1/3$ | Growth exponent: $W(t) \sim t^\beta$ at early times | KPZ (1986) |
| $\alpha$ | $1/2$ | Roughness exponent: $W_{\mathrm{sat}} \sim L^\alpha$ | KPZ (1986) |
| $z$ | $3/2$ | Dynamic exponent: $t_{\times} \sim L^z$ | KPZ (1986) |

Here $W(t) = \sqrt{\langle(h - \langle h\rangle)^2\rangle}$ is the
interface width (roughness).

!!! note "Not equilibrium critical exponents"
    These are **not** the standard $\alpha, \beta$ of thermal
    phase transitions. The KPZ $\beta$ is the growth exponent
    (width vs time), and $\alpha$ is the roughness exponent
    (saturation width vs system size). Do not confuse with
    order-parameter or specific-heat exponents.

### Galilean Invariance Constraint

The non-linear term $(\nabla h)^2$ endows the KPZ equation with
**Galilean invariance** under tilted-frame transformations. This
symmetry enforces the exact relation

$$\alpha + z = 2$$

which, combined with the scaling relation $z = \alpha / \beta$,
fixes all three exponents from a single one:

$$\alpha = 1/2, \quad z = 3/2, \quad \beta = \alpha/z = 1/3$$

### Exact Distribution

Beyond the exponents, the full probability distribution of the
height fluctuations is known exactly in $1+1$D:

- **Flat initial condition**: $\chi \sim t^{1/3} \xi_{\mathrm{GOE}}$
  (Tracy-Widom GOE distribution)
- **Curved initial condition**: $\chi \sim t^{1/3} \xi_{\mathrm{GUE}}$
  (Tracy-Widom GUE distribution)

This was proven rigorously via the connection to the TASEP and
random matrix theory (Sasamoto-Spohn 2010, Amir-Corwin-Quastel
2011).

---

## Higher Dimensions --- Numerical Estimates

For $d \geq 2$ spatial dimensions, no exact solution is known. QAtlas
returns the best-published numerical estimates with their quoted
uncertainties.

### $2+1$D --- Pagnani & Parisi (2015)

| Exponent | Value | Statistical $1\sigma$ | Source |
|----------|-------|------------------------|--------|
| $\beta$ | $0.2415$ | $\pm 0.0015$ | Pagnani–Parisi 2015 |
| $\alpha$ | $0.393$ | $\pm 0.005$ | Pagnani–Parisi 2015 |
| $z$ | $1.613$ | $\pm 0.009$ | Pagnani–Parisi 2015 |

Galilean invariance $\alpha + z = 2$ is satisfied within combined
error bars: $0.393 + 1.613 = 2.006 \pm 0.014$.

### $3+1$D --- Kelling & Ódor (2011)

| Exponent | Value | Estimated $1\sigma$ | Source |
|----------|-------|----------------------|--------|
| $\beta$ | $0.18$ | $\sim 0.01$ | Kelling–Ódor 2011 |
| $\alpha$ | $0.31$ | $\sim 0.01$ | Kelling–Ódor 2011 |
| $z$ | $1.51$ | $\sim 0.01$ | Kelling–Ódor 2011 |

!!! warning "Galilean invariance not strictly satisfied at d=3"
    The numerical estimates above give $\alpha + z \approx 1.82$,
    well below the symmetry-required value of $2.0$. This is a
    *known open issue* in the d≥3 KPZ literature, not a violation
    of Galilean invariance — it reflects the difficulty of
    extracting $\alpha$ and $z$ from finite-size simulations far
    from the d=1 fixed point. Treat the d=3 entry as a best-numerical
    pointer, not as a sharp reference value.

The upper critical dimension of KPZ (if it exists) remains an open
problem; QAtlas does not expose values for $d \geq 4$.

---

## QAtlas API

```julia
using QAtlas

# 1+1D KPZ: exact growth exponents (Rational{Int})
g1 = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=1)
# (β_growth = 1//3, α_rough = 1//2, z = 3//2)

# 2+1D KPZ: numerical estimates with 1σ companions
g2 = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=2)
# (β_growth = 0.2415, β_growth_err = 0.0015,
#  α_rough  = 0.393,  α_rough_err  = 0.005,
#  z        = 1.613,  z_err        = 0.009)

# 3+1D KPZ
g3 = QAtlas.fetch(Universality(:KPZ), GrowthExponents(); d=3)
# (β_growth = 0.18, α_rough = 0.31, z = 1.51, … _err = 0.01 each)

# d ≥ 4 raises an ErrorException.
```

The `_err` fields are absent at $d=1$ (exact) and present at $d=2,3$
(numerical).  Use `haskey(g, :β_growth_err)` to branch on
exact-vs-numerical if needed.

Note the use of `GrowthExponents()` rather than `CriticalExponents()`
to reflect the non-equilibrium nature of the KPZ class.

---

## References

- M. Kardar, G. Parisi, Y.-C. Zhang, "Dynamic scaling of growing
  interfaces", Phys. Rev. Lett. **56**, 889 (1986) --- original KPZ
  equation and exponent prediction.
- M. Prähofer, H. Spohn, "Universal distributions for growth processes
  in 1+1 dimensions and random matrices", Phys. Rev. Lett. **84**,
  4882 (2000) --- d=1 exact distribution.
- T. Sasamoto, H. Spohn, "One-dimensional Kardar-Parisi-Zhang
  equation: an exact solution and its universality",
  Phys. Rev. Lett. **104**, 230602 (2010) --- exact height
  distribution.
- G. Amir, I. Corwin, J. Quastel, "Probability distribution of the
  free energy of the continuum directed random polymer in 1 + 1
  dimensions", Comm. Pure Appl. Math. **64**, 466 (2011).
- I. Corwin, "The Kardar-Parisi-Zhang equation and universality
  class", Random Matrices Theory Appl. **1**, 1130001 (2012) ---
  review.
- A. Pagnani, G. Parisi, "Numerical estimate of the
  Kardar-Parisi-Zhang universality class in (2+1) dimensions",
  Phys. Rev. E **92**, 010101(R) (2015) --- d=2 numerical estimates.
- J. Kelling, G. Ódor, "Extremely large-scale simulation of a
  Kardar-Parisi-Zhang model using graphics cards",
  Phys. Rev. E **84**, 061150 (2011) --- d=3 numerical estimates.
- T. Halpin-Healy, "(2+1)-dimensional directed polymer in a random
  medium: scaling phenomena and universal distributions",
  Phys. Rev. Lett. **109**, 170602 (2012) --- d=2 cross-method
  consistency.
