# Kitaev1D — 1D p-wave Majorana Wire

## Overview

The 1D Kitaev (2001) chain is the canonical exactly-solvable model of a
1D topological superconductor.  It is a free-fermion Bogoliubov-de
Gennes problem with a `Z_2` Pfaffian invariant: the topological phase
hosts two Majorana zero modes at the ends of an open chain.

> **Distinct from `KitaevHoneycomb`.**  This model is the *one-dimensional*
> spinless p-wave wire (Kitaev, *Phys.-Usp.* 44, 131, 2001).  The
> [`KitaevHoneycomb`](kitaev-honeycomb.md) model is a *two-dimensional*
> anisotropic spin model on the honeycomb lattice (Kitaev,
> *Ann. Phys.* 321, 2, 2006).  They share a name only by historical
> accident.

## Hamiltonian

```math
H = -\mu \sum_i c_i^{\dagger} c_i
    - t   \sum_i \bigl(c_i^{\dagger} c_{i+1} + \text{h.c.}\bigr)
    + \Delta \sum_i \bigl(c_i c_{i+1} + \text{h.c.}\bigr)
```

with `c_i` spinless fermions, chemical potential `μ`, hopping `t`, and
p-wave pairing `Δ`.  Defaults: `μ = 0`, `t = 1`, `Δ = 1`.

**PBC dispersion (closed form):**

```math
E(k) = \sqrt{(2t \cos k + \mu)^2 + 4\Delta^2 \sin^2 k}.
```

**Phase diagram** (`Δ ≠ 0`, `t ≠ 0`):

| Regime              | Phase           | Pfaffian invariant |
| ------------------- | --------------- | ------------------ |
| `\|μ\| < 2\|t\|`    | topological     | `ν = -1`           |
| `\|μ\| = 2\|t\|`    | gapless         | ill-defined        |
| `\|μ\| > 2\|t\|`    | trivial         | `ν = +1`           |

**Topological invariant** (Kitaev 2001, Pfaffian sign):

```math
\nu = \operatorname{sgn}\bigl[\operatorname{Pf} A(k=0)
                              \cdot \operatorname{Pf} A(k=\pi)\bigr]
    = \operatorname{sgn}(\mu^2 - 4t^2),
```

evaluated on the 2 × 2 Majorana Bloch matrix `A(k)` at the two
time-reversal-invariant momenta.

**OBC edge zero modes:** in the topological phase the two Majorana ends
hybridise into a single complex fermion with hybridisation energy
`E_edge(N) ~ exp(-N/ξ)` where `ξ ~ 1/log(2|t|/|μ|)` for
`|μ| ≪ 2|t|`.  At the *sweet spot* `μ = 0`, `t = Δ` the two Majoranas
decouple exactly and `E_edge` vanishes for all `N ≥ 2`.

---

## TFIM correspondence

The transverse-field Ising model is a special case of `Kitaev1D` under
the identification

```math
\mu = -2h, \qquad t = J, \qquad \Delta = J.
```

The OBC BdG spectra of `Kitaev1D(μ=-2h, t=J, Δ=J)` and `TFIM(J=J, h=h)`
agree exactly (verified by `test/standalone/test_kitaev1d.jl`).  The
helper `_kitaev1d_bdg_spectrum` is a strict generalisation of
`_tfim_bdg_spectrum`; choosing `Δ = J` and `μ = -2h` reproduces TFIM's
A and B blocks element-wise.

---

## Coverage Matrix

| Quantity                          | OBC                | Infinite                              |
| --------------------------------- | ------------------ | ------------------------------------- |
| [`ExactSpectrum`](@ref)           | BdG (size `N`)     | —                                     |
| [`Energy`](@ref) `{:per_site}`    | conversion         | Gauss-Kronrod over `E(k)`             |
| [`MassGap`](@ref)                 | BdG (smallest)     | analytic min over `k`                 |
| [`EdgeModeEnergy`](@ref)          | BdG (smallest)     | —                                     |
| [`CorrelationLength`](@ref)       | —                  | `1 / Δ_gap` (`Inf` on critical line)  |
| [`TopologicalInvariant`](@ref)    | —                  | Pfaffian sign at `k = 0, π`           |

---

## Usage

```julia
using QAtlas

# Topological phase (μ = 0, sweet spot)
m = Kitaev1D(; μ=0.0, t=1.0, Δ=1.0)
fetch(m, TopologicalInvariant(), Infinite())   # -1
fetch(m, MassGap(), Infinite())                # 2.0  (bulk gap)
fetch(m, EdgeModeEnergy(), OBC(40))            # ~1e-15 (sweet-spot zero modes)

# Trivial phase
m_triv = Kitaev1D(; μ=3.0, t=1.0, Δ=1.0)
fetch(m_triv, TopologicalInvariant(), Infinite())   # +1
fetch(m_triv, MassGap(), Infinite())                # 1.0

# TFIM cross-check
J, h = 1.0, 0.7
m_tfim = Kitaev1D(; μ=-2h, t=J, Δ=J)
fetch(m_tfim, ExactSpectrum(), OBC(20))    # matches _tfim_bdg_spectrum(20, J, h)
```

---

## References

- A. Yu. Kitaev, "Unpaired Majorana fermions in quantum wires",
  *Phys.-Usp.* **44**, 131 (2001).
- J. Alicea, "New directions for the pursuit of Majorana fermions in
  solid state systems", *Rep. Prog. Phys.* **75**, 076501 (2012).
- J. K. Asbóth, L. Oroszlány, A. Pályi,
  *A Short Course on Topological Insulators*,
  Lect. Notes Phys. **919** (2016) — Pfaffian invariant.
