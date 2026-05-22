# Heisenberg Spinons: Dispersion, des Cloizeaux–Pearson Continuum, Müller Ansatz

This page collects the closed-form excitation kinematics of the spin-``\tfrac{1}{2}``
isotropic antiferromagnetic Heisenberg chain in the thermodynamic limit
that are exposed as Phase-1 helpers and `fetch` methods on
[`Heisenberg1D`](@ref).

## Main results

For

```math
H = J\sum_{i}\mathbf{S}_i\cdot\mathbf{S}_{i+1},\qquad J>0,
```

the elementary excitations are massless spinons (half-odd-integer spin)
which always come in pairs in any physical observable. The single-spinon
dispersion (Faddeev–Takhtajan 1981) is

```math
\boxed{\;\varepsilon(k) \;=\; \frac{\pi J}{2}\,|\sin k|,\qquad k\in[0,\pi].\;}
```

The two-spinon continuum, parameterised by the total momentum ``q``, is
bounded by the des Cloizeaux–Pearson (1962) edges

```math
\boxed{\;
  \varepsilon_L(q) \;=\; \frac{\pi J}{2}\,|\sin q|,\qquad
  \varepsilon_U(q) \;=\; \pi J\,\bigl|\sin(q/2)\bigr|.
\;}
```

The lower edge coincides with the single-spinon dispersion and the
continuum is gapless at ``q=0`` and ``q=\pi`` (Umklapp).

The longitudinal dynamic structure factor inside the continuum is
approximated by the Müller ansatz (Müller–Thomas–Beck–Bonner 1981):

```math
\boxed{\;
  S^{zz}_{\rm Müller}(q,\omega)
  \;=\; \frac{\Theta\!\bigl[\omega-\varepsilon_L(q)\bigr]\,
              \Theta\!\bigl[\varepsilon_U(q)-\omega\bigr]}
             {2\,\sqrt{\omega^2 - \varepsilon_L(q)^2}},
\;}
```

with ``S^{zz}=0`` outside ``[\varepsilon_L,\varepsilon_U]``.

---

## Derivation sketch

### Spinon dispersion (Faddeev–Takhtajan 1981)

The Bethe-ansatz solution of the spin-``\tfrac{1}{2}`` XXX antiferromagnet
admits a thermodynamic state — the antiferromagnetic Dirac sea — built
from a continuous distribution of real rapidities ``\lambda\in\mathbb{R}``
with density ``\rho_0(\lambda)=1/(2\cosh(\pi\lambda))``. Holes in this
distribution carry spin ``\tfrac{1}{2}`` — these are the spinons.

The dispersion of a single spinon is obtained by adding one hole of
rapidity ``\lambda`` on top of the sea. The resulting energy and momentum,
relative to the ground state, are

```math
\varepsilon(\lambda) \;=\; \frac{\pi J}{2}\,
                          \frac{1}{\cosh(\pi\lambda)},\qquad
p(\lambda) \;=\; \frac{\pi}{2} \;-\; \arctan\bigl(\sinh(\pi\lambda)\bigr).
```

Eliminating ``\lambda`` via ``\cosh(\pi\lambda) = 1/\sin p`` — itself a
direct consequence of the second relation — yields the closed-form
dispersion stated above:

```math
\varepsilon(p) \;=\; \frac{\pi J}{2}\,\sin p,\qquad p\in[0,\pi].
```

The absolute value ``|\sin p|`` in the boxed formula extends the result by
the periodicity of the Brillouin zone.

### Two-spinon continuum and des Cloizeaux–Pearson edges (1962)

A pair of spinons with momenta ``k_1, k_2 \in [0,\pi]`` carries total
momentum ``q = k_1 + k_2`` (mod ``2\pi``) and total energy
``E = \varepsilon(k_1) + \varepsilon(k_2)``. At fixed ``q`` the
spinon-pair energy ranges over an interval whose endpoints are
extracted by Lagrange-multiplying ``\varepsilon(k_1)+\varepsilon(k_2)``
with the constraint ``k_1 + k_2 = q``:

```math
\partial_{k_1}\varepsilon(k_1) \;=\; \partial_{k_2}\varepsilon(k_2)
  \quad\Longleftrightarrow\quad
  \cos k_1 \;=\; \cos k_2.
```

Two solutions emerge.

* ``k_1 = k_2 = q/2`` — both spinons share the momentum, giving the
  **upper edge**

  ```math
\varepsilon_U(q) \;=\; 2\,\varepsilon(q/2)
                       \;=\; \pi J\,\bigl|\sin(q/2)\bigr|.
```

* ``k_1 = 0,\ k_2 = q`` (one spinon at the gapless point) —
  giving the **lower edge**

  ```math
\varepsilon_L(q) \;=\; \varepsilon(0) + \varepsilon(q)
                       \;=\; \frac{\pi J}{2}\,|\sin q|.
```

Hence ``\varepsilon_L(q) \equiv \varepsilon(q)``, and the continuum
collapses (``\varepsilon_U = \varepsilon_L = 0``) at the gapless points
``q = 0, \pi``.

Numerically, at ``q = \pi`` one has ``\varepsilon_L = 0`` and
``\varepsilon_U = \pi J``, which is the value carried by
[`heisenberg_two_spinon_upper_edge`](@ref) at the AFM ordering wave
vector.

### Müller ansatz for ``S^{zz}(q,\omega)`` (1981)

Müller, Thomas, Beck, and Bonner proposed an explicit closed form for
the longitudinal dynamic structure factor that

* has the correct support on the two-spinon continuum,
* reproduces the integrable square-root singularity at the lower edge
  ``\omega \to \varepsilon_L^+`` (which dominates the spectral weight),
* is normalised to give the correct equal-time longitudinal structure
  factor in leading order.

The ansatz is

```math
S^{zz}_{\rm Müller}(q,\omega) \;=\;
  \frac{\Theta(\omega-\varepsilon_L)\,\Theta(\varepsilon_U-\omega)}
       {2\,\sqrt{\omega^2 - \varepsilon_L^2}},
```

returning ``0`` outside the closed continuum
``[\varepsilon_L(q),\varepsilon_U(q)]``. The ansatz is **approximate**:
it captures the lower-edge behaviour and the support exactly but
misestimates the spectral weight near the upper edge, where four-spinon
contributions become important.

The ansatz value diverges as ``\omega \to \varepsilon_L^+`` but remains
integrable: $\int_{\varepsilon_L}^{\varepsilon_U}
S^{zz}\,d\omega < \infty$. QAtlas returns the raw analytical value
without a regulator; downstream callers integrating in ``\omega`` should
either use a quadrature aware of the square-root singularity (e.g.
Gauss–Chebyshev, or the change of variables ``\omega^2 = \varepsilon_L^2 + s``)
or regulate via ``\sqrt{\omega^2 - \varepsilon_L^2 + \eta^2}`` at their
own choice of ``\eta``.

---

## API

```julia
heisenberg_spinon_dispersion(model::Heisenberg1D, k::Real; J::Real = 1.0)::Float64
heisenberg_two_spinon_lower_edge(model::Heisenberg1D, q::Real; J::Real = 1.0)::Float64
heisenberg_two_spinon_upper_edge(model::Heisenberg1D, q::Real; J::Real = 1.0)::Float64

fetch(::Heisenberg1D, ::ZZStructureFactor, ::Infinite;
      q::Real, ω::Real, method::Symbol = :muller, J::Real = 1.0)::Float64
```

Special values, all at ``J = 1``:

| quantity                              | ``q = 0`` | ``q = \pi/2`` | ``q = \pi`` |
| ------------------------------------- | ------- | ----------- | --------- |
| `heisenberg_spinon_dispersion`        | ``0``     | ``\pi/2``     | ``0``       |
| `heisenberg_two_spinon_lower_edge`    | ``0``     | ``\pi/2``     | ``0``       |
| `heisenberg_two_spinon_upper_edge`    | ``0``     | ``\pi/\sqrt{2}`` | ``\pi``  |

For the Quantity-based dispatch the routing is:

```julia
fetch(Heisenberg1D(), ZZStructureFactor(), Infinite();
      q = π/2, ω = 1.0)                        # → Müller ansatz
fetch(Heisenberg1D(), ZZStructureFactor(), Infinite();
      q = π/2, ω = 1.0, method = :caux_hagemans)
# → Phase-2 placeholder; raises an informative error
```

Quasiparticle dispersion stays a top-level helper (no `QuasiparticleDispersion`
quantity type is introduced in Phase 1) — this matches the existing TFIM
style with [`tfim_quasiparticle_dispersion`](@ref) and
[`tfim_two_spinon_dos`](@ref). A unified `QuasiparticleDispersion`
quantity is a candidate refactor target if/when more models grow
analogous helpers.

---

## Phase 2 (TODO): exact dynamic structure factor

The Müller ansatz is the standard 1981-vintage approximation. The
**exact** result for the dynamic longitudinal structure factor of the
infinite Heisenberg chain is given by the algebraic-Bethe-ansatz form
factor sum due to

> J.-S. Caux, R. Hagemans,
> *The four-spinon dynamical structure factor of the Heisenberg chain*,
> J. Stat. Mech. P12013 (2006).

Implementing the Caux–Hagemans formula is **Phase 2 of issue #154** and
is not yet shipped. The current `fetch(::Heisenberg1D,
::ZZStructureFactor, ::Infinite; method = :caux_hagemans, …)` branch
raises an informative `ErrorException` so downstream code can probe for
availability. The Müller branch (`method = :muller`, default) remains
the only Phase-1 implementation.

---

## References

* J. des Cloizeaux, J. J. Pearson,
  "Spin-wave spectrum of the antiferromagnetic linear chain",
  Phys. Rev. **128**, 2131 (1962).
* L. D. Faddeev, L. A. Takhtajan,
  "What is the spin of a spin wave?",
  Phys. Lett. A **85**, 375 (1981).
* G. Müller, H. Thomas, H. Beck, J. C. Bonner,
  "Quantum spin dynamics of the antiferromagnetic linear chain in zero
  and nonzero magnetic field",
  Phys. Rev. B **24**, 1429 (1981).
* J.-S. Caux, R. Hagemans,
  "The four-spinon dynamical structure factor of the Heisenberg chain",
  J. Stat. Mech. P12013 (2006). (Phase 2 reference.)
