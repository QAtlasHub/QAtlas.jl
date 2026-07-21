# S1Heisenberg1D — Spin-1 AFM Heisenberg Chain (Haldane Chain)

!!! warning "Status: Unstable (v0.18.x)"
    The spin-1 dense-ED observable surface debuted in v0.17. The local
    Hilbert space is 3-dimensional, so the global space is $3^N$ and
    the hard cap is $N \le 8$. Method signatures and kwarg names may
    change in v0.19. `LocalMagnetization(:y)` is intentionally omitted
    (Tier 3 in the observable issue) until a use case appears.

## Hamiltonian

```math
H = J \sum_{i} \mathbf{S}_i \cdot \mathbf{S}_{i+1},
\qquad S = 1\ (\text{3-dimensional local Hilbert space}).
```

The on-site spin operators are the standard $3 \times 3$ matrices
`_S1_x`, `_S1_y`, `_S1_z` with eigenvalues $\{-1, 0, +1\}$.

## Phases

The spin-1 antiferromagnet is in the **Haldane phase**: a gapped,
topologically non-trivial phase with hidden $\mathbb{Z}_2 \times
\mathbb{Z}_2$ string order and (on OBC) free spin-1/2 edge modes.

| Quantity | Value |
|---|---|
| Bulk gap $\Delta_\infty$ | $\approx 0.41048\,J$ (White–Huse 1993, DMRG) |
| GS energy density $e_0$ | $\approx -1.401484\,J$ (White–Huse 1993, DMRG) |
| Edge | gapless spin-1/2 doublet on each end (OBC) |

There is no closed form for $\Delta_\infty$ or $e_0$ — they are exposed
by QAtlas as `:literature_value` rows.

## Coverage Matrix

OBC rows are dense-ED with cap $N \le 8$ ($3^N \le 6561$).

| Quantity | OBC ($N \le 8$) | Infinite |
|---|---|---|
| `Energy{:total}` | dense-ED | — |
| `Energy{:per_site}` | conversion | $\approx -1.40148\,J$ (White–Huse 1993) |
| `FreeEnergy` / `ThermalEntropy` / `SpecificHeat` | dense-ED | — |
| `MagnetizationX` / `Y` / `Z` | dense-ED | — |
| `LocalMagnetization(:x)` / `LocalMagnetization(:z)` | dense-ED | — |
| `SusceptibilityXX` / `YY` / `ZZ` | dense-ED | — |
| `XXCorrelation` / `YY` / `ZZ` (`:static`, `:connected`) | dense-ED | — |
| `VonNeumannEntropy` / `RenyiEntropy` | partial trace | — |
| `MassGap` | dense-ED ($E_1 - E_0$) | $\approx 0.41048\,J$ (Haldane gap) |
| `EnergyLocal` | dense-ED (symmetric bond split) | — |

`LocalMagnetization(:y)` is **not** registered (Tier 3, deferred).

## Spin-1 Convention

The on-site operators are the spin-1 matrices, so eigenvalues run over
$\{-1, 0, +1\}$ (NOT the Pauli $\pm 1$ of TFIM/XXZ1D's $\sigma^\alpha$).
All observables on `S1Heisenberg1D` use the $S^\alpha$ convention
directly. When comparing against TFIM-style outputs that use
$\sigma^\alpha = 2 S^\alpha$, factors of $2$ (linear) or $4$ (variance,
two-point) are required.

## SU(2) Symmetry Identities

`is_su2_symmetric(::S1Heisenberg1D) === true`, so the SU(2) row of
`SYMMETRY_IDENTITIES` (PR #133) is checked automatically:

```math
\chi_{xx} = \chi_{yy} = \chi_{zz},
\qquad m_\alpha = 0\ \ (\alpha \in \{x, y, z\}).
```

This holds to ED precision since the dense kernel never breaks the
spin rotation.

## Code Examples

```julia
using QAtlas

m = S1Heisenberg1D(J=1.0)
N = 4
β = 1.0

QAtlas.fetch(m, Energy(),                     OBC(N); beta=β)
QAtlas.fetch(m, MassGap(),                    OBC(N))             # finite-N gap
QAtlas.fetch(m, MassGap(),                    Infinite())         # ≈ 0.41048
QAtlas.fetch(m, ZZCorrelation{:static}(),     OBC(N); beta=β, i=1, j=4)
QAtlas.fetch(m, VonNeumannEntropy(),          OBC(N); ℓ=2, beta=Inf)
```

For OBC $N \le 4$ the gap $E_1 - E_0$ is dominated by the spin-1/2
edge-mode quasi-degeneracy, so it is much smaller than $\Delta_\infty$.
The two values are not directly comparable until $N$ is large enough
for bulk physics to dominate (typical DMRG convergence around $N \sim
40$).

## References

- F. D. M. Haldane, Phys. Rev. Lett. **50**, 1153 (1983) —
  $\sigma$-model prediction of the integer-spin gap.
- S. R. White, D. A. Huse, Phys. Rev. B **48**, 3844 (1993) —
  state-of-the-art DMRG values for $\Delta_\infty$ and $e_0$.
- T. Kennedy, H. Tasaki, Phys. Rev. B **45**, 304 (1992) — string
  order and hidden $\mathbb{Z}_2 \times \mathbb{Z}_2$ symmetry.
- I. Affleck, J. Phys. Condens. Matter **1**, 3047 (1989) —
  bosonisation of the spin-1 chain.

## Related

- [Heisenberg1D](heisenberg.md) — spin-1/2 cousin; gapless Luttinger
  liquid, topologically trivial — a distinct phase.
- [XXZ1D](xxz.md) — spin-1/2 anisotropic generalisation; SU(2)-symmetric
  only at $\Delta = 1$.
