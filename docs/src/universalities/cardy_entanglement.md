# Calabrese–Cardy Entanglement (Universality-Level API)

## Overview

For any 1+1D conformal field theory of central charge ``c``, the
Calabrese–Cardy formulae give a closed-form universal value for the
von Neumann (and Rényi) entanglement entropy of a contiguous block of
length ``\ell`` on a chain of length ``L``. QAtlas exposes these as
universality-level `fetch` methods so they can be evaluated without
ever instantiating a specific lattice model — the only physics input
is the universality class itself, from which the central charge is
read out via `fetch(Universality{C}(), CentralCharge())`.

This page documents the *universality* dispatch added in
[#149](https://github.com/QAtlasHub/QAtlas.jl/issues/149). The
rigorous twist-operator derivation of the ``c/3`` vs ``c/6`` prefactors
lives separately in
[`docs/src/calc/calabrese-cardy-obc-vs-pbc.md`](../calc/calabrese-cardy-obc-vs-pbc.md).

## Closed forms

Let ``c`` be the central charge of the universality class. Drop the
non-universal additive constant ``c_{1}^{\prime}`` (UV cutoff) and the
Affleck–Ludwig boundary entropy ``\log g`` (boundary-state-dependent),
both of which require model-specific input not available at the
universality level. What remains is the universal log-prefactor
piece:

```math
S_{\rm PBC}(\ell, L) = \frac{c}{3}\,\log\!\left[\frac{L}{\pi}\sin\!\left(\frac{\pi\ell}{L}\right)\right],
```

```math
S_{\rm OBC}(\ell, L) = \frac{c}{6}\,\log\!\left[\frac{2L}{\pi}\sin\!\left(\frac{\pi\ell}{L}\right)\right],
```

```math
S_{\rm \infty}(\ell) = \frac{c}{3}\,\log\ell.
```

The Rényi-``\alpha`` extension uses the substitution

```math
c \;\to\; c\,\frac{1 + 1/\alpha}{2}
```

in the same closed form, which reduces to ``c`` at ``\alpha = 1`` (von
Neumann).

## API

```julia
fetch(::Universality{C}, ::VonNeumannEntropy, ::PBC; ℓ::Real, L::Real, kwargs...)
fetch(::Universality{C}, ::VonNeumannEntropy, ::OBC; ℓ::Real, L::Real, kwargs...)
fetch(::Universality{C}, ::VonNeumannEntropy, ::Infinite; ℓ::Real, kwargs...)

fetch(::Universality{C}, ::RenyiEntropy, ::PBC; ℓ::Real, L::Real, kwargs...)
fetch(::Universality{C}, ::RenyiEntropy, ::OBC; ℓ::Real, L::Real, kwargs...)
fetch(::Universality{C}, ::RenyiEntropy, ::Infinite; ℓ::Real, kwargs...)
```

The central charge is fetched internally via
`fetch(Universality{C}(), CentralCharge(); kwargs...)`. Universality
classes that do not have a 1+1D-CFT central charge defined (e.g. KPZ,
Percolation, the 3D O(``n``) classes) raise `ErrorException` with a
clear message.

## Supported universality classes

| Class | ``d`` | ``c`` | Reference |
|-------|-----|-----|-----------|
| `Universality(:Ising)` | 2 | ``1/2`` | M(3,4) — Belavin–Polyakov–Zamolodchikov, Nucl. Phys. B 241, 333 (1984) |
| `Universality(:IsingSDRG)` | 2 | ``\ln(2)/2`` | Strong-disorder renormalization group (SDRG) / infinite-randomness fixed point (IRFP) of 1D random TFIM — Refael–Moore, Phys. Rev. Lett. 93, 260602 (2004) |
| `Universality(:Potts3)` | 2 | ``4/5`` | M(5,6) — Dotsenko, Nucl. Phys. B 235, 54 (1984) |
| `Universality(:Potts4)` | 2 | ``1`` | Compact boson at marginal point — di Francesco–Mathieu–Sénéchal §12.3 |
| `Universality(:XY)` | 2 | ``1`` | BKT free boson — Kosterlitz J. Phys. C 7, 1046 (1974) |

For all other classes — `Universality(:KPZ)`, `Universality(:Percolation)`,
`Universality(:Heisenberg)` (3D O(3)), and any future class — the
`fetch(::Universality{C}, ::CentralCharge; ...)` method is
intentionally *not* defined, so any Cardy entanglement call routes
through the helper `_cardy_central_charge` and raises
`ErrorException` with the message:

```
Universality{:<C>}: Calabrese-Cardy entanglement requires a 1+1D CFT central
charge, which is not defined for this universality class. Define
fetch(::Universality{:<C>}, ::CentralCharge; ...) first, or use a class that
lives in a 1+1D CFT (e.g. :Ising d=2, :Potts3 d=2, :Potts4 d=2, :XY d=2).
```

## Examples

### Ising critical chain, finite size

```julia
using QAtlas

# c = 1/2; (1/6) log[(8/π) sin(π·4/8)] = (1/6) log(8/π) ≈ 0.1558
S_pbc = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=4.0, L=8.0)
S_obc = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), OBC(); ℓ=4.0, L=8.0)
```

### Thermodynamic limit

```julia
# (1/2)/3 · log(10) = (1/6) log 10
S_inf = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), Infinite(); ℓ=10.0)
```

### Rényi at α = 2

The substitution ``c \to c \cdot (1 + 1/\alpha)/2 = (3/4)c`` at
``\alpha = 2`` makes the Rényi-2 entropy exactly ``3/4`` of the von
Neumann value with the same log argument:

```julia
S_vn = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=4.0, L=8.0)
S_r2 = QAtlas.fetch(Universality(:Ising), RenyiEntropy(2.0), PBC(); ℓ=4.0, L=8.0)
@assert S_r2 ≈ (3//4) * S_vn
```

## Caveats

- **Non-universal constants dropped.** The UV cutoff term ``c_{1}^{\prime}``
  and the Affleck–Ludwig boundary entropy ``\log g`` are *not*
  included in the returned value. Tests that fit lattice
  entanglement-entropy data against this reference must include an
  additive offset.
- **Endpoints.** `ℓ = 0` and `ℓ = L` (the formula limit where the
  block contains zero or all sites) return `-Inf` rather than
  `NaN` or a finite floating-point artefact of `sin(π)`. This is
  the correct continuum limit — the entropy diverges as the cut
  encloses arbitrarily fewer sites and is regularised by the UV
  cutoff that we have dropped.
- **No off-critical extension.** These formulae assume the
  universality class is realised at criticality; gapped universality
  classes do not have a CFT central charge. Off-critical
  entanglement (mass crossover) is the responsibility of
  model-level fetch methods (e.g. TFIM in
  `src/models/quantum/TFIM/TFIM_cft_entanglement.jl`).

## References

- P. Calabrese, J. Cardy, *Entanglement entropy and quantum field
  theory*, J. Stat. Mech. P06002 (2004).
- P. Calabrese, J. Cardy, *Entanglement entropy and conformal field
  theory*, J. Phys. A 42, 504005 (2009).
- I. Affleck, A. W. W. Ludwig, *Universal noninteger "ground-state
  degeneracy" in critical quantum systems*, Phys. Rev. Lett. 67,
  161 (1991) — Affleck–Ludwig ``\log g`` term.
