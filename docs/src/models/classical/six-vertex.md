# Classical Six-Vertex (Ice-Rule) Model on the Square Lattice

## Overview

The six-vertex model is the canonical exactly-solvable lattice model
of *constrained* statistical mechanics. Each edge of the square lattice
carries an arrow; at every vertex the **ice rule** (Pauling 1935;
Lieb 1967a) demands two incoming and two outgoing arrows, leaving
exactly six allowed local configurations. Lieb's three 1967 papers
solved each of the named sub-models in closed form, and Sutherland
1967 followed up with the analytic continuation across the full
disordered phase.

Each of the six vertex configurations is weighted in pairs:

```math
\omega_1 = \omega_2 = a, \qquad \omega_3 = \omega_4 = b, \qquad \omega_5 = \omega_6 = c.
```

All thermodynamic information at fixed temperature is encoded in the
single dimensionless invariant

```math
\Delta = \frac{a^2 + b^2 - c^2}{2 a b}
```

which divides parameter space into three exactly solvable phases:

| Phase                      | Condition       | Status in this version                           |
| -------------------------- | --------------- | ------------------------------------------------ |
| Ferroelectric (FE, frozen) | ``\Delta > 1``    | Closed form (Lieb 1967c)                         |
| Disordered                 | ``|\Delta| \le 1``| Lieb / Sutherland 1967 trigonometric integral    |
| Antiferroelectric (AFE)    | ``\Delta < -1``   | Lieb 1967b elliptic form — *deferred to phase 3* |

The square-ice point ``a = b = c = 1`` (so ``\Delta = 1/2``) sits inside
the disordered phase and admits the celebrated Lieb 1967a residual
entropy.

---

## Square-Ice Residual Entropy (Lieb 1967a)

### Statement

At the symmetric point ``a = b = c = 1`` the zero-temperature
configurational entropy per vertex is the closed form

```math
\frac{S}{N} = \frac{3}{2} \log\frac{4}{3} \approx 0.4315231086776713\ldots
```

### Physical Context

This is the per-vertex residual entropy of two-dimensional
*square ice* — the ground-state ensemble of the ice-rule manifold
without any energy bias. The number of admissible configurations
grows exponentially with ``N`` (the lattice volume) at exactly this rate.
Lieb's derivation uses a Bethe-ansatz solution of the transfer matrix
in the special unit-weight limit; the analogous Pauling 1935 mean-field
estimate ``S/N = \log(3/2)`` is *not* tight.

### API

```julia
m = QAtlas.square_ice()
QAtlas.fetch(m, ResidualEntropy(), Infinite())     # → 0.4315231086776713
```

The implementation returns the closed form `(3/2) * log(4/3)` directly.

---

## Free-Energy Density: Ferroelectric Phase (Lieb 1967c)

### Statement

For ``\Delta > 1`` the ground state is the unique frozen configuration
in which all arrows are parallel along the dominant axis, and

```math
f_{\mathrm{FE}}(a, b, c) = -\log \max(a, b).
```

### Physical Context

The partition function is dominated by a single configuration:
``Z \sim \big(\max(a, b)\big)^N``, so ``f = -\log\max(a, b)``. At the
KDP point ``a > 2``, ``b = c = 1`` this reduces to the Lieb 1967c
result ``f = -\log a``. The phase boundary ``\Delta = 1`` is the KDP
critical point.

### API

```julia
m = QAtlas.kdp_model(2.0)   # a = 2, b = c = 1, Δ > 1 ⇒ FE
QAtlas.fetch(m, FreeEnergy(), Infinite())          # → -log(2)
QAtlas.fetch(m, ResidualEntropy(), Infinite())     # → 0.0  (frozen)
```

---

## Free-Energy Density: Disordered Phase (Lieb / Sutherland 1967)

### Statement

For ``|\Delta| \le 1``, parameterise ``\Delta = -\cos\mu`` with
``0 \le \mu \le \pi``. The free-energy density is the trigonometric
integral

```math
-f(a, b, c) = \log c
            + \frac{1}{\pi} \int_0^\infty
              \frac{\sinh\big((\pi - \mu) x\big)\, \tanh(\mu x)}{x \cosh(\mu x)} \, dx
```

(Lieb 1967a, Sutherland 1967; cf. Baxter 1982 §8.8). The integrand
decays exponentially at large ``x`` and its small-``x`` Taylor expansion
is ``\mu (\pi - \mu) x + \mathcal{O}(x^3)``, so the integral is
well-conditioned.

### Physical Context

Inside the disordered phase the system has algebraic correlations and
a continuously varying critical exponent set as a function of ``\mu``.
The square-ice point is the high-symmetry interior point ``\mu = \pi/3``
(so ``\Delta = 1/2``), and at that special value the integral evaluates
to the closed-form Lieb 1967a result, which we use as a
machine-precision cross-check on the numerical quadrature.

### API

```julia
m = SixVertex(; a=1.0, b=1.0, c=0.5)   # Δ = 0.875, disordered
QAtlas.fetch(m, FreeEnergy(), Infinite())     # finite, evaluated by QuadGK to ~1e-12

# Square-ice cross-check
m = QAtlas.square_ice()
QAtlas.fetch(m, FreeEnergy(), Infinite())     # ≈ -(3/2) log(4/3)
```

The implementation evaluates the integral with `QuadGK.quadgk` on
``(0, \infty)`` with `rtol = 1e-12`, `atol = 1e-14`, using a small-``x``
Taylor expansion for ``|x| < 10^{-8}`` to suppress floating-point
cancellation.

---

## Free-Energy Density: Antiferroelectric Phase (Lieb 1967b) — *Deferred*

For ``\Delta < -1``, parameterise ``\Delta = -\cosh\lambda`` with
``\lambda > 0``. The closed form is the Lieb 1967b
elliptic-function expression involving Jacobi theta functions; it is
heavier than the disordered branch and is **deferred to phase 3** of
issue #163. Calls to `fetch(::SixVertex, ::FreeEnergy, ::Infinite)`
in the AFE phase currently raise an informative `ArgumentError`. The
phase boundary at ``\Delta = -1`` (e.g. f-model with ``c = 2``) is
already covered by the disordered-branch limit ``\mu \to \pi``.

---

## Convenience constructors

```julia
QAtlas.square_ice()          # a = b = c = 1                       (Δ = 1/2)
QAtlas.f_model(c::Real)      # a = b = 1, c free; AFE for c > 2    (Δ = 1 − c²/2)
QAtlas.kdp_model(a::Real)    # b = c = 1, a free; FE for a > 2    (Δ = a/2)
```

These are thin wrappers around `SixVertex(; a, b, c)`.

---

## Verification

The standalone test
[`test/standalone/test_six_vertex.jl`](https://github.com/sotashimozono/QAtlas.jl/blob/main/test/standalone/test_six_vertex.jl)
covers:

- Square ice ``S/N = (3/2) \log(4/3) pprox 0.4315231087`` at
  `atol = 1e-14`.
- Phase classification on ``\Delta = \pm 1`` boundaries
  (KDP at ``a = 2``, F-model at ``c = 2``).
- FE plateau ``f = -\log \max(a, b)`` at the KDP point
  (``a = 3``, ``b = c = 1``) and at a ``b``-dominated point
  (``a = 1``, ``b = 3``, ``c = 1``).
- Square-ice FreeEnergy closed form ``f = -(3/2) \log(4/3)``.
- AFE deferral and generic-disordered deferral both asserted via
  `@test_throws ArgumentError`.

Run it as

```bash
julia --project=test test/standalone/test_six_vertex.jl
```

---

## References

- E. H. Lieb, *Residual entropy of square ice*, **Phys. Rev. 162, 162**
  (1967a). Closed form for the square-ice residual entropy.
- E. H. Lieb, *Exact solution of the F model of an antiferroelectric*,
  **Phys. Rev. Lett. 18, 1046** (1967b). AFE elliptic free energy.
- E. H. Lieb, *Exact solution of the two-dimensional Slater KDP model
  of a ferroelectric*, **Phys. Rev. Lett. 19, 108** (1967c). KDP / FE.
- B. Sutherland, *Exact solution of a two-dimensional model for
  hydrogen-bonded crystals*, **Phys. Rev. Lett. 19, 103** (1967).
  Disordered-phase trigonometric integral.
- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*
  (Academic Press, 1982), ch. 8 — textbook treatment.

---

<!-- ATLAS:HUBS:START -- auto-generated by docs/atlas/generate.jl. Do not edit by hand; edits between these markers are overwritten on next regen. -->

## Verified hubs

In the [Verified Atlas](../../atlas/index.md), this model registers 2 hubs (quantity / BC pair). The badge column shows the R1 assurance level; click a hub link to see the exact `verify(...)` calls, references, and corroboration mechanism.

| Quantity | BC | Assurance | Cards |
|---|---|---|---|
| [`FreeEnergy`](../../atlas/hubs/SixVertex_FreeEnergy_Infinite.md) | `Infinite` | 🟠 uncorroborated-but-feasible | 5 |
| [`ResidualEntropy`](../../atlas/hubs/SixVertex_ResidualEntropy_Infinite.md) | `Infinite` | 🟢 corroborated-at-p | 2 |

<!-- ATLAS:HUBS:END -->














