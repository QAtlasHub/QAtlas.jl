# S1TFIM — Spin-1 Transverse-Field Ising Chain

!!! warning "Status: Unstable (v0.18.x)"
    The spin-1 TFIM dense-ED reference debuted in v0.18 alongside the
    spin-1 XXZ chain (issue #96, Phase 1). The local Hilbert space is
    3-dimensional, so the global space is $3^N$ and the hard cap is
    $N \le 8$. Method signatures and kwarg names may change.

## Hamiltonian

```math
H = -J \sum_{i} S^{z}_{i} S^{z}_{i+1} - h \sum_{i} S^{x}_{i},
\qquad S = 1 \ (\text{3-dimensional local Hilbert space}).
```

The on-site spin operators are the same $3 \times 3$ matrices used by
[`S1Heisenberg1D`](@ref) (eigenvalues $\{-1, 0, +1\}$, $S^{x}$ off-diagonal $1/\sqrt{2}$).
This is **not** the Pauli convention used by the spin-$\tfrac{1}{2}$
[`TFIM`](@ref); cross-comparison with the latter requires factors of $2$
(linear) or $4$ (quadratic).

Unlike the spin-$\tfrac{1}{2}$ TFIM, the spin-1 case is *not* a free-
fermion theory: a single site already lives in a 3-state space and
$(S^{x})^{2}$ is not the identity, so the Jordan–Wigner +
Bogoliubov–de Gennes path used in `TFIM.jl` does not apply. The model
exhibits a non-trivial symmetry-breaking transition at finite $h_c/J$
(numerical, no closed form).

## Coverage Matrix (Phase 1)

OBC rows are dense-ED with cap $N \le 8$ ($3^N \le 6561$). Phase 1 ships
the thermal core matching [`S1Heisenberg1D`](@ref) where applicable:

| Quantity | OBC ($N \le 8$) |
|---|---|
| `Energy{:total}` / `FreeEnergy` / `ThermalEntropy` / `SpecificHeat` | dense-ED |
| `MagnetizationX` / `Y` / `Z` (per-site) | dense-ED |
| `SusceptibilityXX` / `YY` / `ZZ` (per-site) | dense-ED |
| `ZZCorrelation` (`:static`, `:connected`) | dense-ED |
| `MassGap` | dense-ED ($E_1 - E_0$) |
| `ExactSpectrum` | dense-ED (sorted eigenvalues) |

## Phase 2 (deferred)

- $XX$ / $YY$ correlators.
- `MagnetizationXLocal` / `MagnetizationZLocal` (site-resolved).
- `EnergyLocal` (bond-symmetric split).
- `VonNeumannEntropy` / `RenyiEntropy` via partial trace.
- Thermodynamic-limit $h_c/J$ literature constant.

## References

- M. Suzuki, *Prog. Theor. Phys.* **56**, 1454 (1976).
- F. C. Alcaraz and A. L. Malvezzi, *J. Phys. A* **28**, 1521 (1995) — phase diagram.