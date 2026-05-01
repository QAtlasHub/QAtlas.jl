# Heisenberg1D — Spin-1/2 AFM Heisenberg Chain

!!! warning "Status: Unstable (v0.18.x)"
    `Heisenberg1D` is a thin delegator over [XXZ1D](xxz.md) with
    $\Delta = 1$. The struct itself carries **no fields**; the exchange
    coupling $J$ is passed as a kwarg at fetch time. Method signatures
    may change in v0.19. There is also a known wrinkle in the
    `ITensorModels.to_qatlas(::Heisenberg1D)` bridge where $J$ can be
    lost in conversion — for non-unit $J$ prefer `XXZ1D(; J, Δ=1.0)`
    directly.

## Hamiltonian

```math
H = J \sum_{i} \mathbf{S}_i \cdot \mathbf{S}_{i+1},
\qquad \mathbf{S}_i = \tfrac{1}{2}\boldsymbol{\sigma}_i,
\qquad J > 0\ \text{(antiferromagnetic)}.
```

The chain is gapless; the low-energy theory is a $c = 1$ Luttinger
liquid (free compactified boson). The ground state is a singlet
($S_\text{tot} = 0$) for any finite even $N$ with AFM coupling.

## Status: thin delegator over XXZ1D(Δ = 1)

Every OBC observable on `Heisenberg1D` is implemented by forwarding
to `XXZ1D(Δ=1.0, J=J)`:

```julia
QAtlas.fetch(Heisenberg1D(), Energy(), OBC(6); beta=1.0, J=1.5)
# == QAtlas.fetch(XXZ1D(J=1.5, Δ=1.0), Energy(), OBC(6); beta=1.0)
```

For $J = 1$ the kwarg can be omitted (default `1.0`).

## Coverage Matrix

All OBC rows delegate to `XXZ1D(Δ=1)` dense-ED ($N \le 12$). The
infinite-chain ground-state energy density is the original Hulthén
value.

| Quantity | OBC | Infinite |
|---|---|---|
| `Energy` / `FreeEnergy` / `ThermalEntropy` / `SpecificHeat` | via XXZ1D | — |
| `MagnetizationX` / `Y` / `Z` (+ `…Local`) | via XXZ1D | — |
| `SusceptibilityXX` / `YY` / `ZZ` | via XXZ1D | — |
| `XXCorrelation` / `YY` / `ZZ` (`:static`, `:connected`) | via XXZ1D | — |
| `VonNeumannEntropy` / `RenyiEntropy` | via XXZ1D | — |
| `MassGap` | via XXZ1D ED gap | $0$ (gapless Luttinger) |
| `CentralCharge` | — | $1$ (free boson) |
| `GroundStateEnergyDensity` | — | $J(1/4 - \ln 2)$ (Hulthén 1938) |
| `ExactSpectrum` | $N = 2$ OBC dimer + $N = 4$ PBC ring | — |

The closed-form $N = 2$ dimer and $N = 4$ ring spectra are exposed
through `ExactSpectrum` and used as harness anchors.

## SU(2) Symmetry Identities

`is_su2_symmetric(::Heisenberg1D) === true`, so the SU(2) row of
`SYMMETRY_IDENTITIES` (PR #133) is automatically applied by the test
harness. The identities checked numerically are

```math
\chi_{xx} = \chi_{yy} = \chi_{zz},
\qquad m_\alpha = 0\ \ (\alpha \in \{x, y, z\}).
```

These hold to ED precision because the dense-ED kernel diagonalises
the full SU(2)-symmetric Hamiltonian without breaking the rotation.

## Quick-look code

```julia
using QAtlas

m = Heisenberg1D()
β = 1.0
N = 6

QAtlas.fetch(m, Energy(),                OBC(N); beta=β, J=1.0)
QAtlas.fetch(m, SpecificHeat(),          OBC(N); beta=β, J=1.0)
QAtlas.fetch(m, MassGap(),               OBC(N);          J=1.0)
QAtlas.fetch(m, MassGap(),               Infinite();      J=1.0)   # 0
QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite();   J=1.0)   # J(1/4 - ln 2)

# Closed-form spectra (used as test anchors)
QAtlas.fetch(m, ExactSpectrum(); N=2, J=1.0, bc=:OBC)   # dimer
QAtlas.fetch(m, ExactSpectrum(); N=4, J=1.0, bc=:PBC)   # 4-site ring
```

## Closed-form anchors

### Dimer ($N = 2$, OBC)

The two-site Hilbert space splits into a singlet and a triplet:

```math
\text{Spec}(H) = \left\{-\tfrac{3J}{4},\; \tfrac{J}{4},\; \tfrac{J}{4},\; \tfrac{J}{4}\right\}
```

Singlet–triplet gap $\Delta = J$; full derivation in
[Heisenberg dimer: singlet–triplet](../../calc/heisenberg-dimer-singlet-triplet.md).

### 4-site PBC ring

The 16-dimensional spectrum decomposes into a unique singlet ground
state at $E_0 = -2J$, a triplet at $-J$, mixed states at $0$, and a
quintet at $+J$. Used as a finite-size cross-check against the Bethe
ansatz.

### Thermodynamic-limit ground-state energy

```math
e_0 = J\!\left(\tfrac{1}{4} - \ln 2\right) \approx -0.4431\,J
```

Full Bethe-ansatz derivation in
[Bethe ansatz: Heisenberg $e_0$](../../calc/bethe-ansatz-heisenberg-e0.md).

## References

- H. Bethe, Z. Physik **71**, 205 (1931) — original Bethe ansatz.
- L. Hulthén, Ark. Mat. Astron. Fys. **26A**, No. 11 (1938) — $e_0 = J(1/4 - \ln 2)$.
- A. Auerbach, *Interacting Electrons and Quantum Magnetism* (Springer, 1994).
- T. Giamarchi, *Quantum Physics in One Dimension* (Oxford, 2004), Ch. 6.

## Related

- [XXZ1D](xxz.md) — anisotropic generalisation; Heisenberg lives at $\Delta = 1$.
- [S1Heisenberg1D](s1heisenberg.md) — spin-1 cousin in the Haldane phase
  (gapped, topologically non-trivial).
- [TFIM](tfim.md) — $c = 1/2$ Ising critical line for contrast.
