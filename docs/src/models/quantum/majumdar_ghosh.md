# MajumdarGhosh — Spin-½ J₁–J₂ Chain at J₂/J₁ = 1/2

The Majumdar–Ghosh chain is the spin-½ J₁–J₂ Heisenberg chain locked
to the special ratio $J_2/J_1 = 1/2$. At this point the ground state
is *exactly* the product of nearest-neighbour singlets — one of the
cleanest closed-form ground states in 1D quantum magnetism.

## Hamiltonian

```math
H = J \sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+1}
    + \frac{J}{2} \sum_i \mathbf{S}_i \cdot \mathbf{S}_{i+2},
\qquad \mathbf{S}_i = \tfrac{1}{2}\boldsymbol{\sigma}_i,
\qquad J > 0.
```

The next-nearest-neighbour coupling is locked to $J/2$ by the model
definition; only $J$ is a free parameter.

## Exact dimer ground state

At $J_2/J_1 = 1/2$ the ground state is the product of nearest-
neighbour singlets, with two inequivalent dimer coverings (even/odd):

```math
|\psi_0^{\pm}\rangle = \prod_i |\text{singlet}\rangle_{(i, i+1)}.
```

Each singlet contributes $\langle \mathbf{S}\cdot\mathbf{S}\rangle = -3/4$
and adjacent dimers are orthogonal, so all $\mathbf{S}_i\cdot\mathbf{S}_{i+2}$
matrix elements on the dimer state vanish. The size-independent
ground-state energy density is therefore

```math
\boxed{\;\frac{E_0}{N} = -\frac{3J}{8}.\;}
```

The dimer state is an exact eigenstate for both periodic (even $N$)
and open boundary conditions; the ground state is two-fold degenerate
on both.

## Excitation gap

| Source | Value | Method |
|--------|-------|--------|
| Shastry–Sutherland (1981) | $\Delta \geq J/4$ | analytical lower bound |
| White–Affleck (1996) | $\Delta \approx 0.234\,J$ | DMRG (numerical exact) |

Both are exposed through the `MassGap` quantity with a `method` kwarg.

## Coverage Matrix

| Quantity | Infinite | PBC | OBC |
|---|---|---|---|
| `GroundStateEnergyDensity` | $-3J/8$ (analytic) | $-3J/8$ (size-indep.) | — |
| `MassGap` | $J/4$ or $0.234\,J$ | — | — |

## Quick-look code

```julia
using QAtlas

m = MajumdarGhosh(; J=1.0)

QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())            # -3/8
QAtlas.fetch(m, GroundStateEnergyDensity(), PBC(8))                # -3/8
QAtlas.fetch(m, MassGap(), Infinite())                             # 1/4  (Shastry–Sutherland)
QAtlas.fetch(m, MassGap(), Infinite(); method=:lower_bound)        # 1/4
QAtlas.fetch(m, MassGap(), Infinite(); method=:numerical)          # 0.234 (White–Affleck DMRG)
```

## References

- C. K. Majumdar, D. K. Ghosh, "On Next-Nearest-Neighbor Interaction
  in Linear Chain. I/II", *J. Math. Phys.* **10**, 1388 (1969) —
  exact dimer ground state at $J_2/J_1 = 1/2$.
- B. S. Shastry, B. Sutherland, "Excitation spectrum of a dimerized
  next-neighbour antiferromagnetic chain", *J. Phys. C* **14**, L765
  (1981) — analytical lower bound $\Delta \geq J/4$.
- S. R. White, I. Affleck, "Dimerization and incommensurate spiral
  spin correlations in the zigzag spin chain", *Phys. Rev. B* **54**,
  9862 (1996) — DMRG gap $\Delta \approx 0.234\,J$.

## Related

- [Heisenberg1D](heisenberg.md) — the $J_2 = 0$ parent chain (gapless,
  Bethe-ansatz ground state $e_0 = J(1/4 - \ln 2)$).
- [XXZ1D](xxz.md) — anisotropic generalisation of the
  nearest-neighbour Heisenberg point.
