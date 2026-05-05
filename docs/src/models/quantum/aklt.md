# AKLT1D — S=1 Bilinear-Biquadratic Chain at the AKLT Point

!!! warning "Status: Unstable (v0.18.x)"
    `AKLT1D` debuted in v0.18.x as the closed-form "anchor" model for
    the Haldane phase. Method signatures and kwarg names may evolve in
    v0.19. The OBC dense-ED path inherits the spin-1 cap `N ≤ 8` from
    `S1Heisenberg1D`.

## Hamiltonian

```math
H = J \sum_{i} \left[ \mathbf{S}_i \cdot \mathbf{S}_{i+1}
                     + \tfrac{1}{3} (\mathbf{S}_i \cdot \mathbf{S}_{i+1})^2 \right],
\qquad S = 1, \quad J > 0\ \text{antiferromagnetic.}
```

The AKLT point is the special location $\theta = \arctan(1/3)$ on the
$S = 1$ bilinear-biquadratic line where the Hamiltonian becomes — up to
a constant — a sum of bond projectors onto the total-spin-2 subspace:

```math
J \left[ \mathbf{S}_i \cdot \mathbf{S}_{i+1}
        + \tfrac{1}{3}(\mathbf{S}_i \cdot \mathbf{S}_{i+1})^2 \right]
  = 2J\, P_2(i, i+1) - \frac{2J}{3}.
```

Because each bond projector annihilates the Valence Bond Solid (VBS)
state — the matrix product state of bond dimension $D = 2$ built by
splitting each spin-1 site into two spin-$\tfrac{1}{2}$ "halves" and
forming a singlet across every neighbouring pair — the ground state of
$H$ is exactly the VBS, with per-site energy density $e_0 = -2J/3$.

## Phases & Closed-Form Values

`AKLT1D` is the canonical analytical entry point to the **Haldane
phase** of $S = 1$ chains: gapped, topologically non-trivial, with
hidden $\mathbb{Z}_2 \times \mathbb{Z}_2$ symmetry breaking and free
spin-$\tfrac{1}{2}$ edge modes under OBC.

| Quantity | Value | Source |
|---|---|---|
| GS energy density $e_0$ | $-2J/3$ | AKLT 1988 (closed form) |
| Correlation length $\xi$ | $1/\log 3 \approx 0.910$ | AKLT 1988 (closed form) |
| String order parameter $O_{\rm str}$ | $4/9$ | AKLT 1988 + Kennedy-Tasaki 1992 |
| Haldane gap $\Delta_\infty$ | $\approx 0.41048\,J$ | Östlund-Rommer 1995 (DMRG) |
| OBC ground-state degeneracy | 4 (singlet $\oplus$ triplet of edge $\tfrac{1}{2}$-spins) | AKLT 1988 |

The AKLT chain and the spin-1 Heisenberg chain
([`S1Heisenberg1D`](s1heisenberg.md)) are in the **same** topological
phase. They differ in three ways: AKLT has closed-form $e_0$ and
$\xi$; the AKLT bond Hamiltonian is frustration-free under OBC (so the
ground state energy is *exactly* $-(2/3)(N-1)J$ on $N$ sites with no
$1/N^2$ bulk correction); and the gap is essentially the same numerical
value as the spin-1 Heisenberg gap because both models live deep
inside the Haldane phase.

## Coverage Matrix

OBC rows are dense-ED on the $3^N$ Hilbert space, capped at $N \le 8$
($3^8 = 6561$).  All Infinite rows are analytical — closed form for
$e_0$, $\xi$, $O_{\rm str}$, and a literature value for $\Delta$.

| Quantity | Infinite | OBC ($N \le 8$) |
|---|---|---|
| `GroundStateEnergyDensity` | $-2J/3$ (closed form) | — |
| `Energy{:per_site}` | $-2J/3$ (closed form) | conversion via `Energy{:total}` |
| `CorrelationLength` | $1/\log 3$ (closed form) | — |
| `StringOrderParameter` | $4/9$ (closed form) | — |
| `MassGap` | $\approx 0.41048\,J$ (Östlund-Rommer 1995) | — |
| `ExactSpectrum` | — | dense-ED ($3^N \times 3^N$) |

## Usage

```julia
using QAtlas

m = AKLT1D()                         # default J = 1.0
fetch(m, GroundStateEnergyDensity(), Infinite())   # → -2/3
fetch(m, CorrelationLength(),       Infinite())    # → 1/log 3 ≈ 0.910
fetch(m, StringOrderParameter(),    Infinite())    # → 4/9
fetch(m, MassGap(),                 Infinite())    # → 0.41048

# OBC dense ED — full sorted spectrum (3^N entries; N ≤ 8)
λ = fetch(m, ExactSpectrum(), OBC(6))
@assert length(λ) == 3^6
@assert λ[4] - λ[1] < 1e-8                # 4-fold edge degeneracy
@assert λ[1] ≈ -(2/3) * (6 - 1)            # exact frustration-free GS energy
```

## Edge Modes (OBC)

Under open boundary conditions the AKLT chain has two unpaired
spin-$\tfrac{1}{2}$ "halves" — one at each end — that do not get
absorbed into a singlet bond.  These two free $\tfrac{1}{2}$-spins
combine into

```math
\tfrac{1}{2} \otimes \tfrac{1}{2} = 0 \oplus 1
\quad \Longrightarrow \quad
\text{4-fold ground-state degeneracy (singlet} \oplus \text{triplet)}.
```

`fetch(m, ExactSpectrum(), OBC(N))` exhibits this directly: the four
lowest eigenvalues are degenerate to floating-point precision for $N \in
\{4, 6, 8\}$, with the next eigenvalue separated by an $O(J)$ bulk gap
that approaches the Haldane gap as $N$ grows.

## References

- I. Affleck, T. Kennedy, E. H. Lieb, H. Tasaki,
  *"Valence bond ground states in isotropic quantum antiferromagnets"*,
  Commun. Math. Phys. **115**, 477 (1988) — the original AKLT
  construction with $e_0$, $\xi$, $O_{\rm str}$, and the OBC edge
  states.
- I. Affleck, T. Kennedy, E. H. Lieb, H. Tasaki,
  *"Rigorous results on valence-bond ground states in antiferromagnets"*,
  Phys. Rev. Lett. **59**, 799 (1987).
- T. Kennedy and H. Tasaki,
  *"Hidden $\mathbb{Z}_2 \times \mathbb{Z}_2$ symmetry breaking in
  Haldane-gap antiferromagnets"*,
  Phys. Rev. B **45**, 304 (1992) — string order parameter and the
  non-local unitary that exposes it.
- S. Östlund and S. Rommer,
  *"Thermodynamic limit of density matrix renormalization"*,
  Phys. Rev. Lett. **75**, 3537 (1995) — DMRG numerical-exact value of
  the AKLT Haldane gap.
