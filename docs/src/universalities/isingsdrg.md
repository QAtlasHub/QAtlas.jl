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
- **Correlation Lengths (two of them)**: Near criticality ($\delta \propto [\ln J] - [\ln h]$) the *average* correlation length is $\xi \sim \lvert\delta\rvert^{-\nu}$ with $\nu = 2$. The *typical* one is a different, smaller power, $\xi_{\text{typ}} \sim \xi^{1-\psi} \sim \lvert\delta\rvert^{-\nu_{\text{typ}}}$ with $\nu_{\text{typ}} = \nu(1-\psi) = 1$ — the correlation function is not self-averaging, so one length does not describe it.
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

### What this class deliberately does not answer

Every Calabrese–Cardy closed form — the finite-size chord, both Casimir
quantities (`ConformalCasimirEnergy`, `CasimirEnergyCorrection`), Cardy's
density of states, the thermal ``\sinh`` form, the quench light-cone —
**raises an `ErrorException` for `IsingSDRG`**, on every route:

```julia
QAtlas.fetch(Universality(:IsingSDRG), VonNeumannEntropy(), PBC(); ℓ=4.0, L=8.0)
# ERROR: ... the Calabrese-Cardy closed forms are consequences of conformal
#        invariance, and this universality class is not declared to be a 1+1D CFT ...
```

This is a refusal, not a gap in coverage: ``c_{\text{eff}}`` is a logarithmic
coefficient, which is strictly weaker than the conformal invariance those forms
need, and the IRFP scales in an activated way (``\psi = 1/2`` above). Passing
`c` explicitly does not route around it — what is refused is the formula, not
its coefficient.

The ``c_{\text{eff}}/3`` above is the **two-cut** coefficient
(``\text{ncuts}\cdot c_{\text{eff}}/6``), i.e. a block with both edges in the
bulk; one cut takes half.

---

## Exponent table

The exact exponents, and the one quantity this class refuses, generated from the source:

```@autodocs
Modules = [QAtlas]
Pages = ["universalities/IsingSDRG/IsingSDRG.jl"]
Private = false
Order = [:type, :function]
```

---

## References

- D. S. Fisher, "Random transverse field Ising spin chains", Phys. Rev. Lett. **69**, 534 (1992); "Critical behavior of random transverse-field Ising spin chains", Phys. Rev. B **51**, 6411 (1995) --- original SDRG solution and scaling theory.
- G. Refael, J. E. Moore, "Entanglement entropy of random quantum critical points in one dimension", Phys. Rev. Lett. **93**, 260602 (2004) --- derivation of the effective central charge $c_{\text{eff}} = (\ln 2)/2$ at the random-singlet fixed point.
- F. Iglói, C. Monthus, "Strong disorder RG approach of random systems", Phys. Rep. **412**, 277 (2005), [doi:10.1016/j.physrep.2005.02.006](https://doi.org/10.1016/j.physrep.2005.02.006) --- Table 1 (§4.1.2) collects the exponent set returned by `CriticalExponents`; §4.4.2 the Griffiths-phase singularities, §9.1.2 the typical-correlation and Griffiths-exponent relations for general $d$. Equation numbers cited in the source are those of the arXiv version, `cond-mat/0502448`.

---

## Connections

- **Models**: [TFIM](../models/quantum/tfim.md) (with random couplings).
- **Clean counterpart**: [Ising](ising.md) --- clean critical point maps to the $c=1/2$ Ising CFT.
- **Verification**: [Disordered Systems](../verification/disordered.md) --- tests verifying the random-singlet structures.
