# TFIM — Transverse-Field Ising Model

## Overview

The one-dimensional transverse-field Ising model is the canonical
exactly-solvable quantum many-body system with a quantum phase
transition. QAtlas solves it via the
[Jordan-Wigner transformation](../../methods/jordan-wigner/index.md)
followed by Bogoliubov-de Gennes (BdG) diagonalisation; PBC additionally
requires parity-projected (NS+R) sectors.

## Hamiltonian

```math
H = -J \sum_{i} \sigma^z_i \sigma^z_{i+1} - h \sum_{i} \sigma^x_i
```

(σ-convention: eigenvalues ±1.) **Parameters**: Ising coupling $J$
(default 1.0), transverse field $h$.

**Phase diagram**:

- $h/J < 1$: ferromagnetic ordered phase ($\langle \sigma^z\rangle\neq0$)
- $h/J = 1$: quantum critical point ([Ising CFT, $c=1/2$](../../universalities/ising.md))
- $h/J > 1$: quantum paramagnetic phase ($\langle \sigma^z\rangle = 0$)

**Universality**: the critical point belongs to the
[2D Ising universality class](../../universalities/ising.md) via the
quantum-classical mapping (1+1D quantum ↔ 2D classical).

---

## Coverage Matrix

The table below reflects [`TFIM_registry.jl`](https://github.com/Vault/QAtlas.jl/blob/main/src/models/quantum/TFIM/TFIM_registry.jl) `@register` entries.
✅ marks a native fetch method; "conversion" means routed through
another granularity by `core/registry.jl`.

| Quantity                                  | OBC                | PBC                | Infinite                                  |
| ----------------------------------------- | ------------------ | ------------------ | ----------------------------------------- |
| [`Energy`](@ref) `{:total}`               | ✅ BdG             | conversion         | —                                         |
| [`Energy`](@ref) `{:per_site}`            | conversion         | ✅ BdG (NS+R)      | ✅ closed-form                            |
| [`FreeEnergy`](@ref)                      | ✅                 | ✅ NS+R            | ✅                                        |
| [`ThermalEntropy`](@ref)                  | ✅                 | ✅ NS+R            | ✅                                        |
| [`SpecificHeat`](@ref)                    | ✅                 | ✅ NS+R            | ✅                                        |
| [`MagnetizationX`](@ref)                  | ✅                 | ✅ NS+R            | ✅                                        |
| [`MagnetizationZ`](@ref)                  | 0 by Z₂            | —                  | ✅ Pfeuty $m_z = (1-(h/J)^2)^{1/8}$        |
| [`SusceptibilityXX`](@ref)                | ✅ variance        | ✅ NS+R            | ✅ Kubo (Calabrese-Mussardo)              |
| [`SusceptibilityZZ`](@ref)                | ✅ Wick            | —                  | ✅ `N_proxy=80`                            |
| [`CorrelationLength`](@ref)               | —                  | —                  | ✅ $\xi = 1/(2\lvert h-J\rvert)$           |
| [`MassGap`](@ref)                         | ✅                 | ✅                 | ✅ $\Delta = 2\lvert h-J\rvert$            |
| [`XXCorrelation`](@ref) `{:static}`       | ✅ Pfaffian        | —                  | ✅ proxy                                  |
| [`XXCorrelation`](@ref) `{:connected}`    | ✅                 | —                  | ✅                                        |
| [`XXCorrelation`](@ref) `{:dynamic}`      | ✅                 | —                  | —                                         |
| [`ZZCorrelation`](@ref) `{:static}`       | ✅                 | —                  | —                                         |
| [`ZZCorrelation`](@ref) `{:connected}`    | ✅ (= static, Z₂)  | —                  | —                                         |
| [`ZZCorrelation`](@ref) `{:dynamic}`      | ✅                 | —                  | —                                         |
| [`ZZCorrelation`](@ref) `{:lightcone}`    | ✅                 | —                  | —                                         |
| [`ZZStructureFactor`](@ref)               | ✅                 | —                  | ✅ static + dynamic (proxy)               |
| [`VonNeumannEntropy`](@ref)               | ✅ Peschel         | —                  | ✅ CC (T=0 crit/gapped + T>0 crit)         |
| [`RenyiEntropy`](@ref)`(α)`               | ✅ Peschel         | —                  | ✅ CC                                     |
| [`EnergyLocal`](@ref)                     | ✅                 | —                  | —                                         |
| [`MagnetizationXLocal`](@ref) `{:equilibrium}` | ✅                 | —                  | —                                         |
| [`MagnetizationXLocal`](@ref) `{:quench}`     | ✅                 | —                  | ✅ closed-form k-integral (#145)          |
| [`MagnetizationZLocal`](@ref)             | ✅                 | —                  | —                                         |
| [`SpontaneousMagnetization`](@ref)        | —                  | —                  | ✅ alias of `MagnetizationZ`              |
| [`CentralCharge`](@ref)                   | —                  | —                  | ✅ 1/2 (critical) / 0                     |

YY observables (`YYCorrelation`, `SusceptibilityYY`, `MagnetizationY`)
are intentionally **not** implemented — the σʸ JW string makes OBC
contractions expensive; tracked as Tier 3 in issue #110.

---

## Boundary Conditions

QAtlas supports three boundary conditions for the TFIM, each with
different physical content:

| BC       | `fetch` argument | BdG size                   | Physical setting                        |
| -------- | ---------------- | -------------------------- | --------------------------------------- |
| OBC      | `OBC(N)`         | $2N \times 2N$ (numerical) | Open chain, $N$ sites, $N-1$ bonds      |
| PBC      | `PBC(N)`         | parity-projected NS+R      | Ring of $N$ sites, $N$ bonds            |
| Infinite | `Infinite()`     | $k$-integral               | Thermodynamic limit, PBC $N \to \infty$ |

**OBC**: the BdG matrix is diagonalised numerically. Boundary effects
include the Z₂ tunneling splitting in the ordered phase and the
$O(1/N)$ boundary correction at criticality. See
[gap analysis](#energy-gap-and-quantum-phase-transition) below.

**PBC**: the JW transformation produces a fermion parity factor that
splits the partition function into Neveu-Schwarz (anti-periodic) and
Ramond (periodic) sectors with both signs of the parity projector
(LSM). QAtlas evaluates all four (NS±, R±). The Ramond k=0 zero mode
at criticality is handled explicitly.

**Infinite**: the quasiparticle dispersion
$\Lambda(k) = 2\sqrt{J^2 + h^2 - 2Jh\cos k}$ is integrated over the
Brillouin zone using Gauss-Kronrod quadrature (QuadGK.jl).

---

## v0.17 / v0.18 Highlights

!!! warning "Status: Unstable (v0.18.x)"
    The PBC thermodynamics, Z-axis Infinite surface, XX static / connected
    via Pfaffian, Calabrese-Cardy entanglement at Infinite, and dynamic
    structure-factor helpers are **new in v0.17–v0.18**. Method
    signatures, granularity conventions, and keyword-argument names
    (`N_proxy`, `ω`, `ℓ`, `beta`) may change in v0.19. Call sites should
    use the public `QAtlas.fetch(model, quantity, bc; ...)` interface and
    must not depend on internal helpers (the `_tfim_*` prefixed
    functions).

### 1. PBC free-fermion thermodynamics (v0.17)

Jordan-Wigner with a fermion parity factor splits $Z$ into
Neveu-Schwarz and Ramond sectors with both parity-projector signs.
QAtlas sums all four sectors (NS+, NS−, R+, R−); at the critical point
the R-sector $k=0$ zero mode is handled explicitly.

```julia
m  = TFIM(; J=1.0, h=0.5)
β  = 1.0
QAtlas.fetch(m, FreeEnergy(),       PBC(8); beta=β)
QAtlas.fetch(m, MagnetizationX(),   PBC(8); beta=β)
QAtlas.fetch(m, SusceptibilityXX(), PBC(8); beta=β)
QAtlas.fetch(m, MassGap(),          PBC(8))
```

References: Lieb-Schultz-Mattis (1961); Sachdev §4.2.
Source: [`TFIM_pbc_thermal.jl`](https://github.com/Vault/QAtlas.jl/blob/main/src/models/quantum/TFIM/TFIM_pbc_thermal.jl).

### 2. Z-axis Infinite — Pfeuty closed forms (v0.17)

| Quantity                                          | Formula                                               |
| ------------------------------------------------- | ----------------------------------------------------- |
| [`MagnetizationZ`](@ref) (= [`SpontaneousMagnetization`](@ref)) | $m_z = (1 - (h/J)^2)^{1/8}\;\;(h<J)$, else 0           |
| [`CorrelationLength`](@ref)                       | $\xi = 1/(2\lvert h-J\rvert)$ (Inf at criticality)     |
| [`SusceptibilityZZ`](@ref)                        | OBC large-$N$ proxy via `N_proxy` kwarg                |
| [`ZZStructureFactor`](@ref)                       | static $S_{zz}(q)$ from Fourier of large-$N$ correlator|

```julia
QAtlas.fetch(TFIM(; J=1.0, h=0.5), MagnetizationZ(),    Infinite())  # ≈ 0.985
QAtlas.fetch(TFIM(; J=1.0, h=0.7), CorrelationLength(), Infinite())  # 1/0.6 ≈ 1.667
```

Source: [`TFIM_zaxis.jl`](https://github.com/Vault/QAtlas.jl/blob/main/src/models/quantum/TFIM/TFIM_zaxis.jl).

### 3. XX static / connected via Pfaffian Wick (v0.18)

OBC static $\langle\sigma^x_i\sigma^x_j\rangle$ is the $t=0$ limit of
the existing dynamic Wick contraction, evaluated as a real Pfaffian
over the Majorana covariance block. The connected variant subtracts
$\langle\sigma^x_i\rangle\langle\sigma^x_j\rangle$. Infinite uses the
OBC large-$N$ proxy (`N_proxy` kwarg).

```julia
m = TFIM(; J=1.0, h=0.7)
QAtlas.fetch(m, XXCorrelation{:static}(),    OBC(8); beta=Inf, i=3, j=5)
QAtlas.fetch(m, XXCorrelation{:connected}(), OBC(8); beta=Inf, i=3, j=5)
QAtlas.fetch(m, XXCorrelation{:static}(),    Infinite(); i=3, j=5, N_proxy=80)
```

YY OBC remains unimplemented (issue #110, Tier 3).
Source: [`TFIM_xx_static.jl`](https://github.com/Vault/QAtlas.jl/blob/main/src/models/quantum/TFIM/TFIM_xx_static.jl).

### 4. Calabrese-Cardy Infinite entanglement (v0.18)

The thermodynamic-limit von Neumann and Rényi entropies are evaluated
in closed form via the Calabrese-Cardy formula. Coverage:

| Region                   | $T = 0$                                              | $T > 0$                                                              |
| ------------------------ | ---------------------------------------------------- | -------------------------------------------------------------------- |
| Critical (h = J)         | $S = (c/3)\,\log(2\ell)$                             | $S = (c/3)\,\log\!\left[(2\beta/\pi)\sinh(\pi\ell/\beta)\right]$     |
| Gapped (h ≠ J)           | $S = (c/6)\,\log(2\xi\,\sinh(\ell/\xi))$             | error (deferred — see issue #110)                                    |

with $c = 1/2$ for Ising. The Rényi $\alpha\neq 1$ prefactor is
$(c/12)(1 + 1/\alpha)$.

```julia
QAtlas.fetch(TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), Infinite(); ℓ=50)
# ≈ (1/6) log(100) — critical T=0

QAtlas.fetch(TFIM(; J=1.0, h=0.5), RenyiEntropy(2.0),   Infinite(); ℓ=20)
# Rényi-2, gapped CC

QAtlas.fetch(TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), Infinite();
             ℓ=20, beta=4.0)
# critical T>0
```

Source: [`TFIM_cft_entanglement.jl`](https://github.com/Vault/QAtlas.jl/blob/main/src/models/quantum/TFIM/TFIM_cft_entanglement.jl).

### 5. Dynamic structure factor at Infinite (v0.18, proxy)

[`ZZStructureFactor`](@ref) at `Infinite()` is router-dispatched on
the optional `ω` keyword:

- `ω === nothing` → existing static proxy (Fourier of static correlator)
- `ω::Real`       → dynamic proxy (time-evolution + Fourier of dynamic
  correlator)

Two helpers are exported for analytic post-processing:

- [`tfim_quasiparticle_dispersion`](@ref)`(model, k) -> Float64`
  — closed-form Bogoliubov dispersion $\Lambda(k)$.
- [`tfim_two_spinon_dos`](@ref)`(model, ω; q_total = 0.0) -> Float64`
  — two-spinon density of states at fixed total momentum, used to
  identify the continuum threshold.

```julia
m = TFIM(; J=1.0, h=1.0)
QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=π/2, ω=1.5)
tfim_quasiparticle_dispersion(m, π/2)
tfim_two_spinon_dos(m, 1.5; q_total=0.0)
```

Closed-form form-factor expansion (Calabrese-Mussardo) is **not yet
implemented** — issue #110.
Source: [`TFIM_infinite_dynamics.jl`](https://github.com/Vault/QAtlas.jl/blob/main/src/models/quantum/TFIM/TFIM_infinite_dynamics.jl).

---

## Ground-State Energy

### Statement

The ground-state energy of the OBC TFIM with $N$ sites is

```math
E_0 = -\sum_{n=1}^{N} \frac{\Lambda_n}{2}
```

where $\{\Lambda_n\}$ are the positive eigenvalues of the $2N \times 2N$
BdG matrix. At finite temperature $\beta = 1/(k_B T)$:

```math
\langle H \rangle(\beta) = -\sum_{n=1}^{N} \frac{\Lambda_n}{2} \tanh\!\left(\frac{\beta \Lambda_n}{2}\right)
```

### Derivation

The TFIM is solved exactly via the
[Jordan-Wigner transformation](../../methods/jordan-wigner/index.md),
which maps the spin chain to free fermions after a
[Kramers-Wannier duality](../../calc/kramers-wannier-duality.md) step.
The full derivation — including why the duality is needed for the
$\sigma^z\sigma^z$ convention and the explicit construction of the
BdG matrix — is given in the calculation note
**[JW-TFIM-BdG](../../calc/jw-tfim-bdg.md)**.

The result is a $2N \times 2N$ real symmetric BdG matrix whose
eigenvalues come in $\pm\Lambda_n$ pairs. The positive eigenvalues
$\Lambda_n > 0$ are the quasiparticle energies, and the total energy
at inverse temperature $\beta$ is:

```math
\langle H \rangle = -\sum_n \frac{\Lambda_n}{2} \tanh\!\left(\frac{\beta \Lambda_n}{2}\right)
```

!!! note "Thermodynamic limit"
    For PBC in the $N \to \infty$ limit, the quasiparticle dispersion
    is $\Lambda(k) = 2\sqrt{J^2 + h^2 - 2Jh\cos k}$, and the energy
    per site becomes a $k$-integral evaluated by Gauss-Kronrod
    quadrature (QuadGK.jl).

### References

- P. Pfeuty, "The one-dimensional Ising model with a transverse field",
  Ann. Phys. **57**, 79 (1970) — exact solution of the 1D TFIM.
- E. Lieb, T. Schultz, D. Mattis, "Two Soluble Models of an
  Antiferromagnetic Chain", Ann. Phys. **16**, 407 (1961) — JW
  transformation for spin chains.
- S. Sachdev, _Quantum Phase Transitions_, Cambridge University Press
  (2011), Ch. 5 — pedagogical treatment.

### QAtlas API

```julia
m = TFIM(; J=1.0, h=0.5)

# Ground-state energy (β → ∞), OBC, N=16 — total
E₀ = QAtlas.fetch(m, Energy{:total}(), OBC(16))

# Finite-temperature total energy
Eβ = QAtlas.fetch(m, Energy{:total}(), OBC(16); beta=2.0)

# Thermodynamic limit (PBC, N→∞) — per site
ε  = QAtlas.fetch(m, Energy{:per_site}(), Infinite(); beta=2.0)
```

### Verification

| Test file                          | Method                    | What is checked                                        |
| ---------------------------------- | ------------------------- | ------------------------------------------------------ |
| `test_tfim_gap_closure.jl`         | Dense ED via `build_tfim` | $E_0^{\text{ED}} = E_0^{\text{BdG}}$ for $N = 4, 6, 8$ |
| `test_universality_cross_check.jl` | BdG at $N = 200$          | $E_0/N \to -4J/\pi$ at $h = J$                         |

---

## Finite-Temperature Observables

### Statement

At inverse temperature $\beta$ and for $N$ sites (OBC), the following
quantities are computed from the BdG spectrum $\{\Lambda_n\}$:

| Quantity      | Formula                                                                                  | Type                       |
| ------------- | ---------------------------------------------------------------------------------------- | -------------------------- |
| Free energy   | $F = -\frac{1}{\beta}\sum_n \ln\!\left[2\cosh(\beta\Lambda_n/2)\right]$                  | [`FreeEnergy`](@ref)       |
| Entropy       | $S = \beta(\langle H \rangle - F)$                                                       | [`ThermalEntropy`](@ref)   |
| Specific heat | $C_v = -\beta^2\,\partial \langle H \rangle / \partial \beta$                            | [`SpecificHeat`](@ref)     |
| Mag. (X)      | $\langle\sigma^x\rangle$ from the Bogoliubov occupation                                  | [`MagnetizationX`](@ref)   |
| Susc. (XX)    | Variance of $\sum_i \sigma^x_i$ (OBC); Kubo at Infinite                                  | [`SusceptibilityXX`](@ref) |

PBC ⇒ all of the above with parity-projected NS+R sums (v0.17).

### Derivation

All quantities follow from the free-fermion partition function. For
independent modes with energies $\Lambda_n$:

```math
\mathcal{Z} = \prod_n 2\cosh\!\left(\frac{\beta\Lambda_n}{2}\right)
```

The free energy is $F = -\beta^{-1}\ln\mathcal{Z}$, and all other
thermodynamic quantities follow from $\beta$-derivatives.

### References

- S. Sachdev, _Quantum Phase Transitions_ (2011), Ch. 5.3.
- QAtlas: `src/models/quantum/TFIM/TFIM_thermal.jl`,
  `TFIM_pbc_thermal.jl` — full implementation.

### QAtlas API

```julia
m = TFIM(; J=1.0, h=0.5)
β = 2.0

F  = QAtlas.fetch(m, FreeEnergy(),       OBC(16); beta=β)
S  = QAtlas.fetch(m, ThermalEntropy(),   OBC(16); beta=β)
Cv = QAtlas.fetch(m, SpecificHeat(),     OBC(16); beta=β)
Mx = QAtlas.fetch(m, MagnetizationX(),   OBC(16); beta=β)
χ  = QAtlas.fetch(m, SusceptibilityXX(), OBC(16); beta=β)

# PBC variants (NS+R) — v0.17
F_pbc = QAtlas.fetch(m, FreeEnergy(),     PBC(16); beta=β)
Mx_pbc = QAtlas.fetch(m, MagnetizationX(), PBC(16); beta=β)

# Infinite — closed-form k-integrals
F_inf  = QAtlas.fetch(m, FreeEnergy(),     Infinite(); beta=β)
Mx_inf = QAtlas.fetch(m, MagnetizationX(), Infinite(); beta=β)
```

### Verification

| Test file                  | Method                 | What is checked                                       |
| -------------------------- | ---------------------- | ----------------------------------------------------- |
| `test_TFIM_thermal.jl`     | Dense ED ($N \leq 10$) | Exact match of $F$, $S$, $C_v$, $M_x$ vs. ED          |
| `test_TFIM_pbc_thermal.jl` | NS+R vs. ED ($N\leq8$) | PBC parity-projected sums match exact ring partition  |

---

## Energy Gap and Quantum Phase Transition

### Statement

The many-body energy gap $\Delta = E_1 - E_0$ equals the smallest BdG
quasiparticle energy $\Lambda_{\min}$. In the thermodynamic limit:

```math
\Delta = 2|J - h|
```

At the critical point $h = J$, the gap closes as $\Delta \sim N^{-z}$
with dynamic exponent $z = 1$.

### Physical Context

- **Ordered phase** ($h < J$): for OBC with finite $N$, the "gap" seen
  by exact diagonalisation is actually the Z₂ **tunneling splitting**
  between $|\!\uparrow\cdots\uparrow\rangle$ and
  $|\!\downarrow\cdots\downarrow\rangle$, which is exponentially small
  in $N$. This is distinct from the physical excitation gap
  $\Delta \approx 2(J - h)$.
- **Critical point** ($h = J$): $\Delta \sim \pi/N$ (finite-size gap
  for OBC).
- **Disordered phase** ($h > J$): $\Delta \approx 2(h - J)$, the
  paramagnetic gap.

### References

- P. Pfeuty, Ann. Phys. **57**, 79 (1970), Eq. (3.6).
- S. Sachdev, _Quantum Phase Transitions_ (2011), §5.5.

### QAtlas API

```julia
# Infinite chain — closed form Δ = 2|h − J|
QAtlas.fetch(TFIM(; J=1.0, h=0.3), MassGap(), Infinite())   # 1.4
QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), Infinite())   # 0.0  (critical)

# OBC finite-N — smallest positive BdG eigenvalue
QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), OBC(32))      # ≈ π/N

# PBC finite-N — smallest excitation across NS / R sectors (v0.17)
QAtlas.fetch(TFIM(; J=1.0, h=1.0), MassGap(), PBC(16))
```

### Verification

| Test file                          | Method                  | What is checked                                                          |
| ---------------------------------- | ----------------------- | ------------------------------------------------------------------------ |
| `test_tfim_gap_closure.jl`         | Dense ED ($N = 4$–$12$) | Gap shrinks with $N$ at $h = J$                                          |
| `test_tfim_gap_closure.jl`         | ED                      | Ordered-phase gap is Z₂ tunneling ($< 10^{-3}$ for $N=6$)                |
| `test_universality_cross_check.jl` | BdG ($N = 200$)         | $\Delta \approx 2\lvert h-J\rvert$; $\nu z = 1$ from log-log regression  |

---

## Entanglement Entropy at OBC (Peschel)

### Statement

At the critical point $h = J$, the entanglement entropy of a contiguous
block of $\ell$ sites in an $N$-site OBC chain obeys the
[Calabrese-Cardy formula](../../methods/calabrese-cardy/index.md):

```math
S(\ell) = \frac{c}{6}\ln\!\left[\frac{2N}{\pi}\sin\!\left(\frac{\pi \ell}{N}\right)\right] + s_1
```

with central charge $c = 1/2$ (Ising CFT). See the
[Calabrese-Cardy method page](../../methods/calabrese-cardy/index.md)
for OBC vs. PBC prefactors and extraction procedure.

### Physical Context

The TFIM is a free-fermion system after Jordan-Wigner transformation,
so the reduced density matrix on a contiguous block of $\ell$ spins
is Gaussian and its von Neumann (or Rényi) entropy is computable in
$O(\ell^3)$ from the Majorana covariance matrix restricted to that
block (Peschel's correlation-matrix method). QAtlas exposes this
directly via [`VonNeumannEntropy`](@ref) and [`RenyiEntropy`](@ref) at
OBC — no Kramers-Wannier detour is needed, because the internal
$\sigma^x$-string JW convention puts the Majorana pair
$(\gamma_{2i-1}, \gamma_{2i})$ on spin site $i$ directly, and the JW
transformation factorises across any contiguous bipartition up to a
parity factor on $A$ that commutes with $\rho_A$ (Fagotti-Calabrese
2010).

Full derivation of the per-mode entropy formula
$S_A = \sum_k s_2(\nu_k)$ from the Gaussian-preservation theorem, the
Majorana-covariance canonical form, and the contiguous-block JW
factorisation:
**[Peschel correlation-matrix method](../../calc/tfim-entanglement-peschel.md)**.

### References

- P. Calabrese, J. Cardy, J. Stat. Mech. **0406**, P06002 (2004), Eq. (19).
- I. Peschel, J. Phys. A **36**, L205 (2003), Eq. (9).
- G. Vidal, J. I. Latorre, E. Rico, A. Kitaev, Phys. Rev. Lett. **90**,
  227902 (2003).
- M. Fagotti, P. Calabrese, Phys. Rev. Lett. **104**, 227203 (2010).

### QAtlas API

```julia
# Ground-state von Neumann, ℓ = N/2 at criticality
QAtlas.fetch(TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(100); ℓ=50)
# ≈ 0.7256  ((c/6) log((2N/π) sin(πℓ/N)) + s_1, c = 1/2)

# Thermal von Neumann at β = 1
QAtlas.fetch(TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(100); ℓ=50, beta=1.0)

# Rényi α ≠ 1 (v0.18)
QAtlas.fetch(TFIM(; J=1.0, h=1.0), RenyiEntropy(2.0), OBC(100); ℓ=50)
```

### Verification

| Test file                             | Method                 | What is checked                                                                  |
| ------------------------------------- | ---------------------- | -------------------------------------------------------------------------------- |
| `test_TFIM_entanglement.jl`           | Peschel vs. full ED SVD| Machine-precision agreement for every $\ell$ at $N = 10$, three $(J, h)$ points |
| `test_TFIM_entanglement.jl`           | Peschel ($N = 100$)    | Extracted $c \approx 0.5$ within 5% at criticality                              |
| `test_TFIM_entanglement.jl`           | Peschel                | Symmetric $S(\ell) = S(N-\ell)$, area law away from criticality                  |
| `test_TFIM_renyi.jl`                  | Peschel α-trace        | Rényi $\alpha = 2, 3$ matches small-$N$ ED                                       |
| `test_TFIM_cft_entanglement.jl`       | CC at Infinite         | Critical T=0/T>0 and gapped T=0 closed forms vs. analytic                        |
| `test_entanglement_central_charge.jl` | ED ($N \le 14$)        | $c_{\text{extracted}} \approx 0.5$ within 10%                                    |

---

## Coverage by Reference

Physical / methodological backing of each fetch surface:

- **BdG (OBC ground / thermal)**: Pfeuty 1970.
- **PBC parity projection (NS+R)**: Lieb-Schultz-Mattis 1961; Sachdev §4.2.
- **Peschel correlation matrix (entanglement, OBC)**: Peschel 2003;
  Calabrese-Cardy 2004; Fagotti-Calabrese 2010.
- **Calabrese-Cardy (entanglement, Infinite)**: Calabrese-Cardy 2004,
  2009.
- **Pfaffian Wick (XX static / connected)**: Wick 1950 + free-fermion
  Σ contraction.
- **Pfeuty closed forms (Z-axis Infinite)**: Pfeuty 1970 (spontaneous
  magnetisation, correlation length).
- **Two-spinon DOS / dispersion helpers**: standard Bogoliubov
  dispersion + convolution; see Calabrese-Mussardo for the form-factor
  programme (not yet implemented).

---

## Connections

- **Universality**: [Ising universality class](../../universalities/ising.md) —
  $c = 1/2$, $\nu = 1$, $z = 1$.
- **Classical counterpart**: [IsingSquare](../classical/ising-square.md) —
  the 1+1D TFIM maps to the 2D classical Ising model via the
  quantum-classical correspondence
  ($\beta_{\text{classical}} \leftrightarrow$ imaginary time).
- **Disordered version**: [Random TFIM](../../verification/disordered.md) —
  the Fisher infinite-randomness fixed point at
  $[\ln J]_{\text{avg}} = [\ln h]_{\text{avg}}$.
- **E8 spectrum**: [E8 universality](../../universalities/e8.md) —
  perturbing the critical TFIM at $h = J$ by a longitudinal field
  $\lambda \sigma^z$ is the $\Phi_{(1,2)} = \sigma$ magnetic
  perturbation of the Ising CFT. Zamolodchikov (1989) showed the
  resulting massive field theory **remains integrable** and its
  eight stable particles realise the $E_8$ mass spectrum.

## API

`Modules = [QAtlas]` at the end of `index.md` already pulls
docstrings for the exported observable types and TFIM helpers
([`tfim_quasiparticle_dispersion`](@ref),
[`tfim_two_spinon_dos`](@ref)); no `@autodocs` block is needed here.

---

<!-- ATLAS:HUBS:START -- auto-generated by docs/atlas/generate.jl. Do not edit by hand; edits between these markers are overwritten on next regen. -->

## Verified hubs

In the [Verified Atlas](../../atlas/index.md), this model registers 60 hubs (quantity / BC pair). The badge column shows the R1 assurance level; click a hub link to see the exact `verify(...)` calls, references, and corroboration mechanism.

| Quantity | BC | Assurance | Cards |
|---|---|---|---|
| [`CentralCharge`](../../atlas/hubs/TFIM_CentralCharge_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 5 |
| [`ConformalTower`](../../atlas/hubs/TFIM_ConformalTower_OBC.md) | `OBC` | 🟠 uncorroborated-but-feasible | 0 |
| [`ConformalTower`](../../atlas/hubs/TFIM_ConformalTower_PBC.md) | `PBC` | 🟠 uncorroborated-but-feasible | 0 |
| [`CorrelationLength`](../../atlas/hubs/TFIM_CorrelationLength_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 3 |
| [`CriticalExponents`](../../atlas/hubs/TFIM_CriticalExponents_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`Energy`](../../atlas/hubs/TFIM_Energy_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 12 |
| [`Energy`](../../atlas/hubs/TFIM_Energy_OBC.md) | `OBC` | 🟢 corroborated-at-p | 58 |
| [`Energy`](../../atlas/hubs/TFIM_Energy_PBC.md) | `PBC` | 🟢 corroborated-at-p | 16 |
| [`EnergyLocal`](../../atlas/hubs/TFIM_EnergyLocal_OBC.md) | `OBC` | 🟠 uncorroborated-but-feasible | 0 |
| [`FidelitySusceptibility`](../../atlas/hubs/TFIM_FidelitySusceptibility_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 2 |
| [`FidelitySusceptibility`](../../atlas/hubs/TFIM_FidelitySusceptibility_OBC.md) | `OBC` | 🟢 corroborated-at-p | 9 |
| [`FreeEnergy`](../../atlas/hubs/TFIM_FreeEnergy_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 9 |
| [`FreeEnergy`](../../atlas/hubs/TFIM_FreeEnergy_OBC.md) | `OBC` | 🟢 corroborated-at-p | 45 |
| [`FreeEnergy`](../../atlas/hubs/TFIM_FreeEnergy_PBC.md) | `PBC` | 🟢 corroborated-at-p | 30 |
| [`GGEValue`](../../atlas/hubs/TFIM_GGEValue_Infinite.md) | `Infinite` | 🔵 coherent | 1 |
| [`LiebRobinsonBound`](../../atlas/hubs/TFIM_LiebRobinsonBound_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`LiebRobinsonVelocity`](../../atlas/hubs/TFIM_LiebRobinsonVelocity_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`LoschmidtEcho`](../../atlas/hubs/TFIM_LoschmidtEcho_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 6 |
| [`LoschmidtEcho`](../../atlas/hubs/TFIM_LoschmidtEcho_OBC.md) | `OBC` | 🟢 corroborated-at-p | 12 |
| [`MagnetizationX`](../../atlas/hubs/TFIM_MagnetizationX_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 6 |
| [`MagnetizationX`](../../atlas/hubs/TFIM_MagnetizationX_OBC.md) | `OBC` | 🟢 corroborated-at-p | 18 |
| [`MagnetizationX`](../../atlas/hubs/TFIM_MagnetizationX_PBC.md) | `PBC` | 🟢 corroborated-at-p | 18 |
| [`MagnetizationXLocal`](../../atlas/hubs/TFIM_MagnetizationXLocal_Infinite.md) | `Infinite` | 🔵 coherent | 2 |
| [`MagnetizationXLocal`](../../atlas/hubs/TFIM_MagnetizationXLocal_OBC.md) | `OBC` | 🟠 uncorroborated-but-feasible | 0 |
| [`MagnetizationY`](../../atlas/hubs/TFIM_MagnetizationY_OBC.md) | `OBC` | 🟢 corroborated-at-p | 3 |
| [`MagnetizationZ`](../../atlas/hubs/TFIM_MagnetizationZ_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 5 |
| [`MagnetizationZLocal`](../../atlas/hubs/TFIM_MagnetizationZLocal_OBC.md) | `OBC` | 🟠 uncorroborated-but-feasible | 0 |
| [`MassGap`](../../atlas/hubs/TFIM_MassGap_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 25 |
| [`MassGap`](../../atlas/hubs/TFIM_MassGap_OBC.md) | `OBC` | 🟢 corroborated-at-p | 1 |
| [`MassGap`](../../atlas/hubs/TFIM_MassGap_PBC.md) | `PBC` | 🟠 uncorroborated-but-feasible | 0 |
| [`NMRRelaxationExponent`](../../atlas/hubs/TFIM_NMRRelaxationExponent_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`NMRSpinRelaxationRate`](../../atlas/hubs/TFIM_NMRSpinRelaxationRate_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`NMRSpinRelaxationRate`](../../atlas/hubs/TFIM_NMRSpinRelaxationRate_OBC.md) | `OBC` | 🟠 uncorroborated-but-feasible | 0 |
| [`RenyiEntropy`](../../atlas/hubs/TFIM_RenyiEntropy_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`RenyiEntropy`](../../atlas/hubs/TFIM_RenyiEntropy_OBC.md) | `OBC` | 🟢 corroborated-at-p | 74 |
| [`SpecificHeat`](../../atlas/hubs/TFIM_SpecificHeat_Infinite.md) | `Infinite` | 🔵 coherent | 7 |
| [`SpecificHeat`](../../atlas/hubs/TFIM_SpecificHeat_OBC.md) | `OBC` | 🟢 corroborated-at-p | 40 |
| [`SpecificHeat`](../../atlas/hubs/TFIM_SpecificHeat_PBC.md) | `PBC` | 🟢 corroborated-at-p | 18 |
| [`SpontaneousMagnetization`](../../atlas/hubs/TFIM_SpontaneousMagnetization_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 7 |
| [`SusceptibilityXX`](../../atlas/hubs/TFIM_SusceptibilityXX_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`SusceptibilityXX`](../../atlas/hubs/TFIM_SusceptibilityXX_OBC.md) | `OBC` | 🟢 corroborated-at-p | 18 |
| [`SusceptibilityXX`](../../atlas/hubs/TFIM_SusceptibilityXX_PBC.md) | `PBC` | 🟢 corroborated-at-p | 12 |
| [`SusceptibilityYY`](../../atlas/hubs/TFIM_SusceptibilityYY_OBC.md) | `OBC` | 🟢 corroborated-at-p | 18 |
| [`SusceptibilityZZ`](../../atlas/hubs/TFIM_SusceptibilityZZ_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`SusceptibilityZZ`](../../atlas/hubs/TFIM_SusceptibilityZZ_OBC.md) | `OBC` | 🟢 corroborated-at-p | 18 |
| [`ThermalEntropy`](../../atlas/hubs/TFIM_ThermalEntropy_Infinite.md) | `Infinite` | 🔵 coherent | 11 |
| [`ThermalEntropy`](../../atlas/hubs/TFIM_ThermalEntropy_OBC.md) | `OBC` | 🟢 corroborated-at-p | 41 |
| [`ThermalEntropy`](../../atlas/hubs/TFIM_ThermalEntropy_PBC.md) | `PBC` | 🟢 corroborated-at-p | 15 |
| [`VonNeumannEntropy`](../../atlas/hubs/TFIM_VonNeumannEntropy_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 0 |
| [`VonNeumannEntropy`](../../atlas/hubs/TFIM_VonNeumannEntropy_OBC.md) | `OBC` | 🟢 corroborated-at-p | 42 |
| [`XXCorrelation`](../../atlas/hubs/TFIM_XXCorrelation_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 21 |
| [`XXCorrelation`](../../atlas/hubs/TFIM_XXCorrelation_OBC.md) | `OBC` | 🟢 corroborated-at-p | 1 |
| [`XXStructureFactor`](../../atlas/hubs/TFIM_XXStructureFactor_Infinite.md) | `Infinite` | 🔵 coherent | 14 |
| [`XXStructureFactor`](../../atlas/hubs/TFIM_XXStructureFactor_OBC.md) | `OBC` | 🟢 corroborated-at-p | 27 |
| [`YYCorrelation`](../../atlas/hubs/TFIM_YYCorrelation_OBC.md) | `OBC` | 🟢 corroborated-at-p | 1 |
| [`YYStructureFactor`](../../atlas/hubs/TFIM_YYStructureFactor_Infinite.md) | `Infinite` | 🔵 coherent | 14 |
| [`YYStructureFactor`](../../atlas/hubs/TFIM_YYStructureFactor_OBC.md) | `OBC` | 🟢 corroborated-at-p | 27 |
| [`ZZCorrelation`](../../atlas/hubs/TFIM_ZZCorrelation_OBC.md) | `OBC` | 🟢 corroborated-at-p | 3 |
| [`ZZStructureFactor`](../../atlas/hubs/TFIM_ZZStructureFactor_Infinite.md) | `Infinite` | 🔵 coherent | 14 |
| [`ZZStructureFactor`](../../atlas/hubs/TFIM_ZZStructureFactor_OBC.md) | `OBC` | 🔵 coherent | 24 |

<!-- ATLAS:HUBS:END -->
























---

<!-- ATLAS:DOCS:START -- auto-generated by docs/atlas/generate.jl. Do not edit by hand; edits between these markers are overwritten on next regen. -->

## API

Every `fetch(::Model, …)` method registered for this model — together with the model struct(s) and exported helpers — generated directly from the source (in lock-step with `@register`):

```@autodocs
Modules = [QAtlas]
Pages = ["models/quantum/TFIM/TFIM.jl", "models/quantum/TFIM/TFIM_cft_entanglement.jl", "models/quantum/TFIM/TFIM_dynamics.jl", "models/quantum/TFIM/TFIM_entanglement.jl", "models/quantum/TFIM/TFIM_fidelity.jl", "models/quantum/TFIM/TFIM_gge.jl", "models/quantum/TFIM/TFIM_infinite_dynamics.jl", "models/quantum/TFIM/TFIM_local.jl", "models/quantum/TFIM/TFIM_loschmidt.jl", "models/quantum/TFIM/TFIM_pbc_thermal.jl", "models/quantum/TFIM/TFIM_quench_entanglement.jl", "models/quantum/TFIM/TFIM_registry.jl", "models/quantum/TFIM/TFIM_sigma_x_quench.jl", "models/quantum/TFIM/TFIM_thermal.jl", "models/quantum/TFIM/TFIM_xx_static.jl", "models/quantum/TFIM/TFIM_xx_yy_structure_factor.jl", "models/quantum/TFIM/TFIM_yy.jl", "models/quantum/TFIM/TFIM_zaxis.jl"]
Private = false
Order = [:type, :function]
```

<!-- ATLAS:DOCS:END -->
