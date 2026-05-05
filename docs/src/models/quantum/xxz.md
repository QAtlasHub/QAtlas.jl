# XXZ1D — Spin-1/2 XXZ Chain

!!! warning "Status: Unstable (v0.18.x)"
    The XXZ1D OBC dense-ED full observable surface was introduced in
    v0.17. Method signatures and kwarg names (`beta`, `i`, `j`, `ℓ`, …)
    may change in v0.19. The infinite-system finite-temperature surface
    is currently restricted to the **XX point (Δ = 0)** via the
    free-fermion (Jordan-Wigner) integrals; general-Δ TBA / NLIE is
    tracked in issue #108.

## Hamiltonian

```math
H = J \sum_{i} \bigl[\, S^x_i S^x_{i+1} + S^y_i S^y_{i+1}
                       + \Delta\, S^z_i S^z_{i+1} \,\bigr]
```

with $\mathbf{S}_i = \tfrac{1}{2}\boldsymbol{\sigma}_i$, exchange
coupling $J$ (default `1.0`), and anisotropy $\Delta$ (default `0.0`,
the XX point).

## Phases

| Regime | Phase | Closed form on the GS energy |
|--------|-------|------------------------------|
| $\Delta < -1$ | Gapped ferromagnet (Ising-like FM) | — (TBA, deferred) |
| $\Delta = -1$ | Saturated ferromagnet | $e_0/J = -1/4$ |
| $-1 < \Delta < 1$ | Luttinger liquid, $c = 1$ | — (general $\Delta$ TBA) |
| $\Delta = 0$ | XX / free fermion | $e_0/J = -1/\pi$ |
| $\Delta = 1$ | Isotropic AF Heisenberg | $e_0/J = 1/4 - \ln 2$ (Hulthén 1938) |
| $\Delta > 1$ | Gapped Néel AFM | — (TBA, deferred) |

## Coverage Matrix

OBC rows are dense-ED (Hilbert dim $2^N$, cap $N \le 12$). Infinite
rows are analytic / Bethe-ansatz closed forms.

| Quantity | OBC | Infinite |
|---|---|---|
| `Energy{:total}` (any $\Delta$) | dense-ED | — |
| `Energy{:per_site}` | conversion via $E/N$ | GS at $\Delta \in \{-1, 0, 1\}$; finite-T at Δ = 0 (issue #108 for general Δ) |
| `FreeEnergy` / `ThermalEntropy` / `SpecificHeat` | dense-ED | Δ = 0 free-fermion (QuadGK); other Δ NaN+warn (issue #108) |
| `MagnetizationX` / `Y` / `Z` (+ `…Local`) | dense-ED | — |
| `SusceptibilityXX` / `YY` / `ZZ` | variance | — |
| `XXCorrelation` / `YY` / `ZZ` (`:static`, `:connected`) | dense-ED | — |
| `VonNeumannEntropy` / `RenyiEntropy` | partial trace | — |
| `MassGap` | dense-ED ($E_1 - E_0$) | $0$ on $-1 < \Delta \le 1$, `NaN` otherwise |
| `CentralCharge` | — | $1$ on critical regime, `NaN` otherwise |
| `LuttingerParameter` | — | $K = \pi / (2(\pi - \gamma))$, $\gamma = \arccos\Delta$ |
| `LuttingerVelocity` | — | $u = (\pi J / 2)\,\sin\gamma / \gamma$ |
| `EnergyLocal` | dense-ED (symmetric bond split) | — |

`SpinWaveVelocity` is a type-level alias of `LuttingerVelocity`.

## XX Point (Δ = 0) — Free-Fermion Thermo at Infinite()

After Jordan-Wigner the XX chain is non-interacting; the
single-particle dispersion (in the spin convention `Sᵅ = σᵅ/2`) is

```math
\varepsilon(k) = -J \cos k, \quad k \in [-\pi, \pi].
```

QAtlas exposes per-site `FreeEnergy`, `Energy{:per_site}`,
`ThermalEntropy`, and `SpecificHeat` at `Infinite()` for `Δ = 0` via
adaptive Gauss-Kronrod quadrature on `[0, π]`:

```math
\begin{aligned}
f(\beta) &= -\frac{1}{\pi\beta} \int_0^\pi \log\!\left(2\cosh\tfrac{\beta\varepsilon(k)}{2}\right) dk, \
e(\beta) &= -\frac{1}{2\pi} \int_0^\pi \varepsilon(k)\,\tanh\tfrac{\beta\varepsilon(k)}{2}\,dk, \
s(\beta) &= \beta\,\bigl(e(\beta) - f(\beta)\bigr), \
C(\beta) &=  \frac{1}{\pi} \int_0^\pi \bigl(\tfrac{\beta\varepsilon}{2}\bigr)^2 \operatorname{sech}^2\!\tfrac{\beta\varepsilon}{2}\,dk.
\end{aligned}
```

```julia
m = XXZ1D(J=1.0, Δ=0.0)
QAtlas.fetch(m, Energy(),         Infinite(); beta=10.0)   # → ≈ -1/π = -0.3183
QAtlas.fetch(m, FreeEnergy(),     Infinite(); beta=1.0)
QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=0.01)   # → ≈ log 2
QAtlas.fetch(m, SpecificHeat(),   Infinite(); beta=1.0)    # → > 0
```

For any Δ ≠ 0 these four calls emit a `@warn` and return `NaN` —
the general-Δ thermal Bethe ansatz / NLIE is tracked in
[issue #108](https://github.com/sotashimozono/QAtlas.jl/issues/108).

## v0.17 Highlights — Dense-ED Full Suite

A single `_xxz1d_thermal_kernel(model, N, β)` performs one
eigendecomposition of the $2^N \times 2^N$ Hamiltonian and reuses the
spectrum / eigenvectors across every observable on the OBC row. The
hard cap is `N ≤ 12` (Hilbert dimension `2^12 = 4096`).

```julia
using QAtlas

m = XXZ1D(J=1.0, Δ=0.5)
β = 1.0

QAtlas.fetch(m, FreeEnergy(),         OBC(6); beta=β)
QAtlas.fetch(m, SpecificHeat(),       OBC(6); beta=β)
QAtlas.fetch(m, MagnetizationZ(),     OBC(6); beta=β)         # = 0  (U(1) conservation)
QAtlas.fetch(m, ZZCorrelation{:static}(), OBC(6); beta=β, i=2, j=4)
QAtlas.fetch(m, RenyiEntropy(2.0),    OBC(6); ℓ=3, beta=Inf)
```

### Δ = 1 isotropic point: SU(2) symmetry identities

At $\Delta = 1$ the chain is SU(2)-symmetric, so every observable
satisfies the harness-checked identities

```math
\chi_{xx} = \chi_{yy} = \chi_{zz}, \qquad m_\alpha = 0\ \ (\alpha \in \{x,y,z\}).
```

```julia
m_iso = XXZ1D(J=1.0, Δ=1.0)
QAtlas.fetch(m_iso, SusceptibilityXX(), OBC(6); beta=1.0)  # =
QAtlas.fetch(m_iso, SusceptibilityYY(), OBC(6); beta=1.0)  # =
QAtlas.fetch(m_iso, SusceptibilityZZ(), OBC(6); beta=1.0)
```

These three calls return the same value to ED precision, and are
checked by the SU(2) row of `SYMMETRY_IDENTITIES` (PR #133).

## Critical-regime CFT data

For $-1 < \Delta < 1$ the chain flows to a $c = 1$ compactified-boson
CFT with continuously varying compactification radius. QAtlas exposes
the standard triplet:

```julia
QAtlas.fetch(XXZ1D(; Δ=0.3), CentralCharge(),     Infinite())  # → 1.0
QAtlas.fetch(XXZ1D(; Δ=0.0), LuttingerParameter(), Infinite())  # → 1.0
QAtlas.fetch(XXZ1D(; Δ=1.0), LuttingerVelocity(),  Infinite())  # → π/2
QAtlas.fetch(XXZ1D(; Δ=1.5), CentralCharge(),     Infinite())  # → NaN (+ warn)
```

Full derivation: [XXZ Luttinger parameters from Bethe
ansatz](../../calc/xxz-luttinger-parameters.md).

## References

- H. Bethe, Z. Physik **71**, 205 (1931).
- L. Hulthén, Ark. Mat. Astron. Fys. **26A**, No. 11 (1938) — $\Delta = 1$ value.
- C. N. Yang, C. P. Yang, Phys. Rev. **150**, 321 (1966) — general-$\Delta$
  Bethe-ansatz integral equation.
- M. Takahashi, *Thermodynamics of One-Dimensional Solvable Models*
  (Cambridge UP, 1999), Ch. 4.
- T. Giamarchi, *Quantum Physics in One Dimension* (Oxford, 2004), Ch. 6.

## Related

- [Heisenberg1D](heisenberg.md) — isotropic SU(2) point ($\Delta = 1$); a
  thin delegator to `XXZ1D(Δ=1.0)`.
- [S1Heisenberg1D](s1heisenberg.md) — spin-1 generalisation; gapped
  Haldane phase, distinct universality.
- [TFIM](tfim.md) — $c = 1/2$ Ising chain for contrast with the
  $c = 1$ XXZ critical line.
