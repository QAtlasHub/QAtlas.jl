# CFT Finite-Size Casimir Correction

## Overview

For a 1+1D conformal field theory at a critical point with central
charge $c$ and CFT velocity $v$, the ground-state energy on a finite
system of size $L$ admits a universal $1/L$ correction
(Cardy 1986; Blöte–Cardy–Nightingale 1986; Affleck 1986):

$$E_0^{\mathrm{PBC}}(L) = L\,\varepsilon_\infty - \frac{\pi c v}{6 L} + O(L^{-2})$$

$$E_0^{\mathrm{OBC}}(L) = L\,\varepsilon_\infty + \varepsilon_{\mathrm{surf}} - \frac{\pi c v}{24 L} + O(L^{-2})$$

The PBC term comes from quantising the CFT on a cylinder of
circumference $L$; the OBC term comes from the corresponding strip.
Their ratio is exactly **4**, a kinematic consequence of the conformal
map.

QAtlas exposes only the universal $1/L$ correction term — the
[`CasimirEnergyCorrection`](@ref) quantity — via
`Universality{C}` dispatch.  The extensive piece
$L\,\varepsilon_\infty$ and the OBC surface term
$\varepsilon_{\mathrm{surf}}$ are model-specific and live on the model
side (e.g. TFIM ground-state energy at the critical field).

---

## API

```julia
using QAtlas

# 1+1D Ising at criticality, periodic chain of size L = 16, v = 2J = 2
QAtlas.fetch(Universality(:Ising), CasimirEnergyCorrection(), PBC(); L=16.0, v=2.0)
# -> -π/96  ≈ -0.0327249...

# Same Ising, OBC
QAtlas.fetch(Universality(:Ising), CasimirEnergyCorrection(), OBC(); L=16.0, v=2.0)
# -> -π/384 ≈ -0.0081812...

# Heisenberg chain (SU(2)_1 WZW; v = (π/2) J at the AFM point)
QAtlas.fetch(Universality(:Heisenberg), CasimirEnergyCorrection(), PBC(); L=16.0, v=π/2)
```

The CFT velocity `v` is **model-dependent** and supplied by the caller:

| Model                       | Velocity $v$ at criticality                |
|----------------------------|--------------------------------------------|
| TFIM, $h = J$              | $v = 2J$                                   |
| AFM Heisenberg chain        | $v = (\pi/2)\,J$                           |
| XXZ Luttinger liquid        | $v = v_F = \pi J \sin(\gamma)/\gamma$ (Bethe ansatz)|

QAtlas already exposes `LuttingerVelocity` / `FermiVelocity` /
`SpinWaveVelocity` for the relevant models — read those, then pass
the value through as the `v` kwarg.

---

## Supported Universality Classes

| `C`            | $c$       | 1+1D CFT origin                                  |
|----------------|-----------|--------------------------------------------------|
| `:Ising`       | $1/2$     | Virasoro minimal model $\mathcal{M}(3,4)$        |
| `:Potts3`      | $4/5$     | Virasoro minimal model $\mathcal{M}(5,6)$        |
| `:Potts4`      | $1$       | Free-boson radius limit (marginal)               |
| `:XY`          | $1$       | Compact free boson / 1+1D Luttinger liquid       |
| `:Heisenberg`  | $1$       | $\mathrm{SU}(2)_1$ Wess–Zumino–Witten model      |

Other classes (`:KPZ`, `:Percolation`, `:MeanField`, …) raise
`ErrorException`:

- **KPZ** is non-equilibrium and has no CFT central charge in the
  Cardy sense.
- **Percolation** has $c = 0$ but the underlying CFT is *logarithmic*
  (non-unitary); the simple Cardy formula does not apply in the same
  form.
- **Mean-field** lives above the upper critical dimension, with no
  1+1D CFT representation.

---

## Verification properties

- **PBC : OBC ratio = 4** (kinematic, class-independent):

  $$\frac{E_{0,\mathrm{Casimir}}^{\mathrm{PBC}}}{E_{0,\mathrm{Casimir}}^{\mathrm{OBC}}} = \frac{-\pi c v / (6 L)}{-\pi c v / (24 L)} = 4.$$

- **$L \to \infty$** $\Rightarrow$ value $\to 0$ as $1/L$.

- **Sign**: always negative (the conformal cylinder lowers the
  ground-state energy below the bulk value).

These properties are exercised in
`test/standalone/test_universality_cft_casimir.jl`.

---

## Phase 2 (TODO)

The conformal **tower of states** —

$$E_n - E_0 = \frac{2\pi v}{L}\,(h_n + \bar h_n) + O(L^{-2})$$

with primary scaling dimensions $(h, \bar h)$ — is tracked separately
as future work (issue #150 Phase 2) and will be exposed via a
`ConformalTower` quantity once implemented.  For 2D Ising the primary
spectrum is

$$\{(0,0),\ (1/16, 1/16),\ (1/2, 1/2)\}$$

(identity, spin field $\sigma$, energy field $\varepsilon$); see
[Ising universality](ising.md) for the corresponding scaling
dimensions.

---

## References

1. J. Cardy, "Operator content of two-dimensional conformally
   invariant theories", *Nucl. Phys. B* **270**, 186 (1986).
2. H. W. J. Blöte, J. L. Cardy, M. P. Nightingale, "Conformal
   invariance, the central charge, and universal finite-size
   amplitudes at criticality", *Phys. Rev. Lett.* **56**, 742 (1986).
3. I. Affleck, "Universal term in the free energy at a critical point
   and the conformal anomaly", *Phys. Rev. Lett.* **56**, 746 (1986).
