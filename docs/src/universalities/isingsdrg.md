# Ising SDRG Universality Class

## Overview

The `IsingSDRG` universality class describes the critical behavior of the one-dimensional random transverse-field Ising model (RTFIM) at its infinite-randomness fixed point (IRFP). Rather than being characterized by a standard conformal field theory (CFT), the scaling behavior is governed by the **Strong-Disorder Renormalization Group (SDRG)** (or Ma-Dasgupta-Hu-Fisher procedure), where the ground state is an ensemble of random singlets (or random spin-pairs).

**Symmetry**: $\mathbb{Z}_2$ (spin-flip symmetry).

**Models in this class**: 1D random transverse-field Ising chain at critical tuning $[\ln J] = [\ln h]$.

**Entanglement Entropy**: The entanglement entropy scales logarithmically, but the prefactor is determined by the **Refael-Moore effective central charge** $c_{\text{eff}} = \ln(2)/2 \approx 0.34657$ rather than the clean CFT value $c = 1/2$.

---

## Properties at the Infinite-Randomness Fixed Point (IRFP)

At the critical point of the random chain, the properties differ drastically from the clean Ising CFT:

- **Activated Dynamic Scaling**: The characteristic energy scale $\Omega$ scales with length scale $L$ as $\ln \Omega \sim L^\psi$, with the critical exponent $\psi = 1/2$.
- **Effective Central Charge**: The disorder-averaged entanglement entropy of a block of size $\ell$ scales as:
  $$\overline{S(\ell)} = \frac{c_{\text{eff}}}{3} \ln \ell + \text{const}$$
  with $c_{\text{eff}} = \frac{\ln 2}{2} \approx 0.34657359$.
- **Typical Correlation Length**: Near criticality ($\delta \propto [\ln J] - [\ln h]$), the correlation length scales as $\ln \xi \sim \lvert\delta\rvert^{-\nu}$ with exponent $\nu = 2$.
- **Fluctuations**: The entanglement entropy has non-vanishing sample-to-sample variance even in the thermodynamic limit.

---

## QAtlas API

In QAtlas, the effective central charge of this universality class can be queried at the universality level:

```julia
using QAtlas

# Query Refael-Moore effective central charge
c_eff = QAtlas.fetch(Universality(:IsingSDRG), CentralCharge(); d=2)
# => 0.34657359027997264 (log(2.0)/2.0)
```

---

## References

- D. S. Fisher, "Random transverse field Ising spin chains", Phys. Rev. Lett. **69**, 534 (1992); "Critical behavior of random transverse-field Ising spin chains", Phys. Rev. B **51**, 6411 (1995) --- original SDRG solution and scaling theory.
- G. Refael, J. E. Moore, "Entanglement entropy of random quantum critical points in one dimension", Phys. Rev. Lett. **93**, 260602 (2004) --- derivation of the effective central charge $c_{\text{eff}} = (\ln 2)/2$ at the random-singlet fixed point.

---

## Connections

- **Models**: [TFIM](../models/quantum/tfim.md) (with random couplings).
- **Clean counterpart**: [Ising](ising.md) --- clean critical point maps to the $c=1/2$ Ising CFT.
- **Verification**: [Disordered Systems](../verification/disordered.md) --- tests verifying the random-singlet structures.
