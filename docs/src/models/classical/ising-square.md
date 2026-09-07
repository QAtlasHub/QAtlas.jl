# Classical 2D Ising Model on the Square Lattice

## Overview

The two-dimensional Ising model on the square lattice is arguably the
most studied model in statistical physics. It was the first model to
exhibit a rigorous second-order phase transition (Onsager, 1944) and
remains the benchmark for critical phenomena.

$$H = -J \sum_{\langle i,j \rangle} \sigma_i \sigma_j, \qquad \sigma_i \in \{-1, +1\}$$

**Parameters**: Ising coupling $J > 0$ (ferromagnetic).

**Universality**: The critical point belongs to the
[2D Ising universality class](../../universalities/ising.md) with
central charge $c = 1/2$.

---

## Partition Function (Transfer Matrix)

### Statement

For an $L_x \times L_y$ square lattice with periodic boundary conditions
in both directions:

$$Z(L_x, L_y, \beta, J) = \mathrm{Tr}(T^{L_x})$$

where $T$ is the $2^{L_y} \times 2^{L_y}$ symmetric transfer matrix.
The matrix elements are

$$T_{\sigma, \sigma'} = e^{\beta J E_h(\sigma)/2} \cdot e^{\beta J E_v(\sigma, \sigma')} \cdot e^{\beta J E_h(\sigma')/2}$$

with $E_h(\sigma) = \sum_{j=1}^{L_y} \sigma_j \sigma_{(j \bmod L_y)+1}$
(horizontal bonds within a row, PBC) and
$E_v(\sigma, \sigma') = \sum_{j=1}^{L_y} \sigma_j \sigma'_j$ (vertical
bonds between rows).

### Physical Context

- Valid for any finite $L_x \times L_y$ with PBC
- $\beta = 0$: $Z = 2^{L_x L_y}$ (all configurations equally weighted)
- $J = 0$: $Z = 2^{L_x L_y}$ (no interactions)

### Derivation

The transfer matrix formulation rewrites the partition function as a
product of row-to-row Boltzmann weights. Each "transfer" from row
$\sigma$ to row $\sigma'$ includes the vertical bonds between the two
rows and half of the horizontal bonds within each row (symmetric split).

The computation uses $Z = \mathrm{tr}(T^{L_x})$ via matrix
exponentiation rather than eigendecomposition, to support
[automatic differentiation](../../methods/automatic-differentiation/index.md)
with `ForwardDiff.Dual` numbers.

Full derivation of the symmetric transfer matrix from the
row-by-row factorisation:
**[Transfer matrix: symmetric split construction
](../../calc/transfer-matrix-symmetric-split.md)**.

Full derivation of the AD path from $\ln Z$ to thermodynamic
observables ($F$, $\langle E\rangle$, $C_v$, $S$, …) via
forward-mode dual numbers:
**[Thermodynamic quantities from $\partial \ln Z$ via AD
](../../calc/ad-thermodynamics-from-z.md)**.

### References

- L. Onsager, "Crystal Statistics. I. A Two-Dimensional Model with an
  Order-Disorder Transition", Phys. Rev. **65**, 117 (1944).
- B. M. McCoy and T. T. Wu, *The Two-Dimensional Ising Model*,
  Harvard University Press (1973), Ch. 2.

### QAtlas API

```julia
Z = QAtlas.fetch(IsingSquare(), PartitionFunction();
                 Lx=4, Ly=4, β=0.44, J=1.0)
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_ising_2x2_classical.jl` | Brute-force $2^N$ enumeration | $Z_{\text{TM}} \approx Z_{\text{BF}}$ for $2\times2$, $2\times3$, $3\times3$ |
| `test_ising_ad_thermodynamics.jl` | ForwardDiff | $\langle E \rangle = -\partial(\ln Z)/\partial\beta$, $C_v = \beta^2 \partial^2(\ln Z)/\partial\beta^2$ |
| `test_ising_square_pfaffian.jl` | Special values | $\beta = 0 \Rightarrow Z = 2^N$, monotonicity, $L_x \leftrightarrow L_y$ symmetry |

---

## Critical Temperature (Onsager)

### Statement

$$T_c = \frac{2J}{\ln(1 + \sqrt{2})} \approx 2.269\,J$$

Equivalently, the critical reduced coupling $K_c = J / T_c$ satisfies

$$\sinh(2K_c) = 1$$

which is the self-dual point of the Kramers-Wannier duality.

### Derivation

The Kramers-Wannier duality maps the high-temperature expansion of $Z$
to the low-temperature expansion of the dual partition function. The
critical point is the unique fixed point of this duality:

$$e^{-2K_c} = \tanh K_c \quad \Longleftrightarrow \quad \sinh(2K_c) = 1$$

Solving: $K_c = \frac{1}{2}\mathrm{arcsinh}(1) = \frac{1}{2}\ln(1 + \sqrt{2})$.

Full derivation of the Kramers–Wannier duality (classical + TFIM
operator form) and the self-dual fixed point:
**[Kramers–Wannier duality
](../../calc/kramers-wannier-duality.md)**.

### References

- H. A. Kramers, G. H. Wannier, "Statistics of the Two-Dimensional
  Ferromagnet. Part I", Phys. Rev. **60**, 252 (1941) — duality.
- L. Onsager, Phys. Rev. **65**, 117 (1944) — exact solution confirming
  the duality prediction.

### QAtlas API

```julia
Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature(); J=1.0)
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_ising_onsager_yang.jl` | Numerical value | $T_c \approx 2.2692$ |
| `test_ising_onsager_yang.jl` | Duality identity | $\sinh(2\beta_c J) = 1$ |
| `test_universality_cross_check.jl` | Yang $M(T_c) = 0$ | Phase boundary consistency |

---

## Spontaneous Magnetization (Yang)

### Statement

$$M(T) = \begin{cases}
\left(1 - \sinh^{-4}(2\beta J)\right)^{1/8} & T < T_c \\
0 & T \geq T_c
\end{cases}$$

The exponent $1/8$ directly gives the order parameter critical exponent
$\beta = 1/8$ of the [Ising universality class](../../universalities/ising.md).

### Derivation

Yang's calculation (1952) uses the Pfaffian method to evaluate the
spontaneous magnetization of the infinite square lattice. The key step
is computing $\langle \sigma_0 \rangle$ as a Toeplitz determinant, which
in the thermodynamic limit gives the closed-form expression above.

Full derivation of the Toeplitz-determinant reduction, the strong
Szegő-theorem evaluation, and the emergence of the exponent
$1/8$ from the Wiener–Hopf $1/2$ × Szegő $1/4$ × correlator
$1/2$ factor chain:
**[Yang's spontaneous magnetization via Toeplitz determinant
](../../calc/yang-magnetization-toeplitz.md)**.

### Physical Context

- $T = 0$ ($\beta \to \infty$): $M = 1$ (fully ordered)
- $T \to T_c^-$: $M \sim (T_c - T)^{1/8}$ (critical exponent $\beta = 1/8$)
- $T \geq T_c$: $M = 0$ (disordered)

!!! note "Critical exponent β = 1/8"
    The exponent $1/8$ is very small, meaning the magnetization
    vanishes slowly: $M \approx 0.7$ even at 1% below $T_c$.

### References

- C. N. Yang, "The spontaneous magnetization of a two-dimensional Ising
  model", Phys. Rev. **85**, 808 (1952).

### QAtlas API

```julia
M = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization();
                 β=0.5, J=1.0)
```

### Verification

| Test file | Method | What is checked |
|-----------|--------|-----------------|
| `test_ising_onsager_yang.jl` | Special values | $M(\beta \to \infty) = 1$, $M(T_c) = 0$ |
| `test_ising_onsager_yang.jl` | Monotonicity | $M$ increases as $T \to 0$ |
| `test_ising_onsager_yang.jl` | $\beta$ extraction | log-log slope $\to 1/8$ near $T_c$ |
| `test_universality_cross_check.jl` | Cross-check | Extracted $\beta$ matches `Universality(:Ising).β` |

---

## Connections

- **Universality**: [Ising universality class](../../universalities/ising.md) —
  $\beta = 1/8$, $\nu = 1$, $c = 1/2$.
- **Quantum counterpart**: [TFIM](../quantum/tfim.md) — the 2D classical
  Ising model maps to the 1+1D TFIM via quantum-classical correspondence.
- **Transfer matrix method**: [Methods](../../methods/transfer-matrix/index.md) —
  computational details of the $\mathrm{tr}(T^{L_x})$ evaluation.
- **AD verification**: [Automatic Differentiation](../../methods/automatic-differentiation/index.md) —
  thermodynamic quantities from $\partial(\ln Z)/\partial\beta$.

### Underlying derivations

Every stored value on this page is backed by a step-by-step
derivation note under `docs/src/calc/`:

- **Partition function** → [`transfer-matrix-symmetric-split`](../../calc/transfer-matrix-symmetric-split.md)
  (row-by-row factorisation, symmetric-split construction of $T$,
  $Z = \mathrm{Tr}(T^{L_x})$).
- **Thermodynamic quantities from $Z$** → [`ad-thermodynamics-from-z`](../../calc/ad-thermodynamics-from-z.md)
  (derivatives of $\ln Z$ via ForwardDiff, fluctuation-dissipation
  identity).
- **Critical temperature** → [`kramers-wannier-duality`](../../calc/kramers-wannier-duality.md)
  (self-dual fixed point $\sinh(2K_c) = 1$).
- **Spontaneous magnetisation** → [`yang-magnetization-toeplitz`](../../calc/yang-magnetization-toeplitz.md)
  (Toeplitz determinant + strong Szegő theorem → $\beta = 1/8$).
- **Critical exponents from CFT** → [`ising-cft-primary-operators`](../../calc/ising-cft-primary-operators.md)
  (Kac table of $\mathcal{M}(3, 4)$; three primaries and their
  scaling dimensions).
- **Scaling relations among exponents** → [`ising-scaling-relations`](../../calc/ising-scaling-relations.md)
  (Rushbrooke, Widom, Fisher, Josephson — all four from the Widom
  scaling form).

---

<!-- ATLAS:HUBS:START -- auto-generated by docs/atlas/generate.jl. Do not edit by hand; edits between these markers are overwritten on next regen. -->

## Verified hubs

In the [Verified Atlas](../../atlas/index.md), this model registers 13 hubs (quantity / BC pair). The badge column shows the R1 assurance level; click a hub link to see the exact `verify(...)` calls, references, and corroboration mechanism.

| Quantity | BC | Assurance | Cards |
|---|---|---|---|
| [`CriticalExponents`](../../atlas/hubs/IsingSquare_CriticalExponents_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`CriticalTemperature`](../../atlas/hubs/IsingSquare_CriticalTemperature_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 6 |
| [`Energy`](../../atlas/hubs/IsingSquare_Energy_Infinite.md) | `Infinite` | 🔵 coherent | 1 |
| [`Energy`](../../atlas/hubs/IsingSquare_Energy_PBC.md) | `PBC` | 🟢 corroborated-at-p | 4 |
| [`FreeEnergy`](../../atlas/hubs/IsingSquare_FreeEnergy_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`FreeEnergy`](../../atlas/hubs/IsingSquare_FreeEnergy_PBC.md) | `PBC` | 🟢 corroborated-at-p | 2 |
| [`PartitionFunction`](../../atlas/hubs/IsingSquare_PartitionFunction_PBC.md) | `PBC` | 🟢 corroborated-at-p | 24 |
| [`SpecificHeat`](../../atlas/hubs/IsingSquare_SpecificHeat_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`SpecificHeat`](../../atlas/hubs/IsingSquare_SpecificHeat_PBC.md) | `PBC` | 🟢 corroborated-at-p | 4 |
| [`SpontaneousMagnetization`](../../atlas/hubs/IsingSquare_SpontaneousMagnetization_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 4 |
| [`ThermalEntropy`](../../atlas/hubs/IsingSquare_ThermalEntropy_Infinite.md) | `Infinite` | 🔵 coherent | 1 |
| [`ThermalEntropy`](../../atlas/hubs/IsingSquare_ThermalEntropy_PBC.md) | `PBC` | 🟢 corroborated-at-p | 4 |
| [`UniversalityClass`](../../atlas/hubs/IsingSquare_UniversalityClass_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |

<!-- ATLAS:HUBS:END -->




























---

<!-- ATLAS:DOCS:START -- auto-generated by docs/atlas/generate.jl. Do not edit by hand; edits between these markers are overwritten on next regen. -->

## API

Every `fetch(::Model, …)` method registered for this model — together with the model struct(s) and exported helpers — generated directly from the source (in lock-step with `@register`):

```@autodocs
Modules = [QAtlas]
Pages = ["models/classical/IsingSquare/IsingSquare.jl", "models/classical/IsingSquare/IsingSquare_registry.jl", "models/classical/IsingSquare/IsingSquare_thermal.jl"]
Private = false
Order = [:type, :function]
```

<!-- ATLAS:DOCS:END -->
