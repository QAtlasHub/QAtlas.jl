# Classical 2D Ising Model on the Triangular Lattice

## Overview

The classical Ising model on the triangular lattice is the canonical
example of *frustrated* statistical mechanics. With antiferromagnetic
coupling each elementary triangle cannot simultaneously satisfy all
three bonds, leading to a macroscopically degenerate ground-state
manifold and a non-zero zero-temperature entropy per site (Wannier
1950). With ferromagnetic coupling the lattice supports a standard
order-disorder transition with ``T_c = 4|J|/\ln 3`` (Houtappel 1950).

```math
H = +J \sum_{\langle i,j \rangle} \sigma_i \sigma_j, \qquad \sigma_i \in \{-1, +1\}
```

(Wannier 1950 sign convention; each site has six nearest neighbours.)

**Parameters**: Ising coupling ``J``. ``J > 0`` — antiferromagnetic
(frustrated). ``J < 0`` — ferromagnetic.

---

## Critical Temperature

### Statement

Wannier 1950 / Houtappel 1950 closed forms:

```math
T_c =
\begin{cases}
0, & J > 0 \quad \text{(AFM, frustrated)} \\
\dfrac{4 |J|}{\ln 3} \approx 3.6410\,|J|, & J < 0 \quad \text{(FM)}
\end{cases}
```

### Physical context

- **AFM (``J > 0``)**: The triangular plaquette is the prototypical
  frustrated unit. The ground-state manifold has extensive entropy
  (see [Residual Entropy](#residual-entropy-wannier-1950) below) and
  no spontaneous symmetry breaking occurs at any ``T > 0``. Spin
  correlations decay algebraically (Stephenson 1964) — the system is
  effectively critical for all ``T > 0``.
- **FM (``J < 0``)**: Standard 2D Ising universality class with
  exponents ``\beta = 1/8``, ``\nu = 1``, ``\eta = 1/4`` (same as
  [`IsingSquare`](ising-square.md)).

### References

- G. H. Wannier, "Antiferromagnetism. The triangular Ising net",
  Phys. Rev. **79**, 357 (1950).
- R. M. F. Houtappel, "Order-disorder in hexagonal lattices",
  Physica **16**, 425 (1950).

### QAtlas API

```julia
# AFM (frustrated): T_c = 0
Tc_afm = QAtlas.fetch(IsingTriangular(; J=1.0), CriticalTemperature(), Infinite())

# FM (Houtappel): T_c = 4 |J| / log 3
Tc_fm = QAtlas.fetch(IsingTriangular(; J=-1.0), CriticalTemperature(), Infinite())
```

### Verification

| Test file | What is checked |
| --- | --- |
| `test_ising_triangular.jl` | AFM `T_c = 0`; FM `T_c = 4|J|/log 3` to 1e-12; ``|J|``-scaling |

---

## Residual Entropy (Wannier 1950)

### Statement

For the antiferromagnetic case (``J > 0``), Wannier (1950) showed by
exact transfer-matrix evaluation that the zero-temperature entropy per
site of the triangular Ising net equals

```math
\frac{S}{N k_B} = \frac{2}{\pi} \int_0^{\pi/3} \ln(2 \cos\theta)\, d\theta \approx 0.32306594722\ldots
```

This is strictly between ``0`` and ``\ln 2 \approx 0.693`` — frustration
admits exponentially many ground states, but not all ``2^N``
configurations.

For the ferromagnetic case (``J < 0``) the ground-state manifold consists
only of the two ferromagnetically polarised states related by the
global ``\mathbb{Z}_2`` spin flip, so ``S_\text{residual} = 0`` in the
thermodynamic limit.

### Physical context

- The Wannier integral evaluates to
  ``S/N \approx 0.32306594722``, consistent with a residual ground-state
  degeneracy of ``\Omega \approx (e^{0.3231})^N \approx 1.381^N`` — a
  finite fraction ``\approx 0.4663`` of ``\log 2``.
- The integrand ``\ln(2 \cos\theta)`` is smooth on ``[0, \pi/3]`` with
  ``2 \cos(\pi/3) = 1``, so the upper endpoint is regular and `QuadGK`
  reaches machine-precision quadrature.

### References

- G. H. Wannier, Phys. Rev. **79**, 357 (1950).
- R. M. F. Houtappel, Physica **16**, 425 (1950) — independent
  derivation in the same period; the kagome-lattice closed form is
  obtained by the same method.

### QAtlas API

```julia
S = QAtlas.fetch(IsingTriangular(; J=1.0), ResidualEntropy(), Infinite())
# 0.32306594722...
```

`fetch` evaluates the Wannier integral via `QuadGK.quadgk` with
`rtol = atol = 1e-14`, so the returned value is accurate to roughly
``10^{-12}``.

### Verification

| Test file | What is checked |
| --- | --- |
| `test_ising_triangular.jl` | ``S/N`` matches ``0.32306594722`` at 1e-9 and the QuadGK recomputation at 1e-12 |
| `test_ising_triangular.jl` | ``J``-independence of ``S/N`` in the AFM branch |
| `test_ising_triangular.jl` | ``0 < S/N < \ln 2`` |

---

## Future work

- **Two-point correlations** ``\langle \sigma_0 \sigma_R \rangle``:
  Stephenson (J. Math. Phys. **5**, 1009, 1964) gave the exact
  asymptotic forms (algebraic decay along the symmetry axes for the
  AFM, exponential for the FM). Tracked as a follow-up issue.
- **Free energy density at finite ``T``**: Wannier 1950 / Houtappel 1950
  give the exact integral representation; not yet wired into a
  `FreeEnergy` fetch method.
- **Kagome-lattice analogue**: same Houtappel 1950 method gives the
  closed form for kagome-Ising; tracked separately.

---

## Connections

- [`IsingSquare`](ising-square.md) — the non-frustrated square-lattice
  counterpart (Onsager 1944, Yang 1952).
- [Ising universality class](../../universalities/ising.md) — relevant
  for the FM branch of `IsingTriangular`. The AFM branch is
  effectively critical for all ``T > 0`` and does *not* sit at a single
  RG fixed point with these exponents.

---

<!-- ATLAS:HUBS:START -- auto-generated by docs/atlas/generate.jl. Do not edit by hand; edits between these markers are overwritten on next regen. -->

## Verified hubs

In the [Verified Atlas](../../atlas/index.md), this model registers 9 hubs (quantity / BC pair). The badge column shows the R1 assurance level; click a hub link to see the exact `verify(...)` calls, references, and corroboration mechanism.

| Quantity | BC | Assurance | Cards |
|---|---|---|---|
| [`CriticalExponents`](../../atlas/hubs/IsingTriangular_CriticalExponents_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`CriticalTemperature`](../../atlas/hubs/IsingTriangular_CriticalTemperature_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 6 |
| [`Energy`](../../atlas/hubs/IsingTriangular_Energy_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`FreeEnergy`](../../atlas/hubs/IsingTriangular_FreeEnergy_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`ResidualEntropy`](../../atlas/hubs/IsingTriangular_ResidualEntropy_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 2 |
| [`SpecificHeat`](../../atlas/hubs/IsingTriangular_SpecificHeat_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`SpontaneousMagnetization`](../../atlas/hubs/IsingTriangular_SpontaneousMagnetization_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`ThermalEntropy`](../../atlas/hubs/IsingTriangular_ThermalEntropy_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`ZZCorrelation`](../../atlas/hubs/IsingTriangular_ZZCorrelation_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |

<!-- ATLAS:HUBS:END -->























---

<!-- ATLAS:DOCS:START -- auto-generated by docs/atlas/generate.jl. Do not edit by hand; edits between these markers are overwritten on next regen. -->

## API

Every `fetch(::Model, …)` method registered for this model — together with the model struct(s) and exported helpers — generated directly from the source (in lock-step with `@register`):

```@autodocs
Modules = [QAtlas]
Pages = ["models/classical/IsingTriangular/IsingTriangular.jl", "models/classical/IsingTriangular/IsingTriangular_registry.jl", "models/classical/IsingTriangular/IsingTriangular_thermal.jl"]
Private = false
Order = [:type, :function]
```

<!-- ATLAS:DOCS:END -->
