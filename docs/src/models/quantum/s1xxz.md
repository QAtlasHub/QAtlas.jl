# S1XXZ1D — Spin-1 XXZ Chain

!!! warning "Status: Unstable (v0.18.x)"
    The spin-1 XXZ dense-ED reference debuted in v0.18 alongside
    [`S1TFIM`](@ref) (issue #96, Phase 1). The local Hilbert space is
    3-dimensional, so the global space is $3^N$ and the hard cap is
    $N \le 8$. Method signatures and kwarg names may change.

## Hamiltonian

```math
H = J \sum_{i} \left[ S^{x}_{i} S^{x}_{i+1} + S^{y}_{i} S^{y}_{i+1}
                    + \Delta\, S^{z}_{i} S^{z}_{i+1} \right],
\qquad S = 1.
```

At $\Delta = 1$ this coincides with [`S1Heisenberg1D`](@ref) (Haldane chain);
at $\Delta = 0$ it is the spin-1 XY model. For large positive $\Delta$ the
chain has Néel order; for $\Delta \lesssim -1$ it enters a ferromagnetic
phase (Schulz 1986). Phase 1 exposes finite-$N$ exact dense-ED only — no
thermodynamic-limit closed forms are claimed here.

The spin-1 operators carry eigenvalues $\{-1, 0, +1\}$, **distinct** from
the Pauli convention of the spin-$\tfrac{1}{2}$ [`XXZ1D`](@ref).

## Coverage Matrix (Phase 1)

OBC rows are dense-ED with cap $N \le 8$ ($3^N \le 6561$):

| Quantity | OBC ($N \le 8$) |
|---|---|
| `Energy{:total}` / `FreeEnergy` / `ThermalEntropy` / `SpecificHeat` | dense-ED |
| `MagnetizationX` / `Y` / `Z` (per-site) | dense-ED |
| `SusceptibilityXX` / `YY` / `ZZ` (per-site) | dense-ED |
| `XXCorrelation` / `YY` / `ZZ` (`:static`, `:connected`) | dense-ED |
| `MassGap` | dense-ED ($E_1 - E_0$) |
| `ExactSpectrum` | dense-ED (sorted eigenvalues) |

The spin-1 Heisenberg literature constants ($e_0 \approx -1.40148\,J$,
$\Delta_\infty \approx 0.41048\,J$) are exposed under [`S1Heisenberg1D`](@ref)
since `S1XXZ1D(; Δ = 1)` numerically reproduces the same Hamiltonian.

## Phase 2 (deferred)

- Phase-diagram literature constants for $\Delta \ne 1$.
- Local observables / partial-trace entropies.
- Thermodynamic-Bethe-ansatz spin-1 free energy at general $\Delta$.

## References

- H. J. Schulz, *Phys. Rev. B* **34**, 6372 (1986) — phase diagram.
- F. D. M. Haldane, *Phys. Lett. A* **93**, 464 (1983).