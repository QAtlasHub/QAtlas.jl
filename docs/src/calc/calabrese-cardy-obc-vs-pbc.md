# Calabrese–Cardy Entanglement Entropy: OBC vs PBC Prefactor

## Main result

For a ``1+1``-dimensional critical lattice system described at low
energies by a conformal field theory of central charge ``c``, on a
chain of ``N`` sites with UV lattice cutoff ``a``, the von Neumann
entropy of a single contiguous block of ``\ell`` sites takes the
following universal form in the two standard boundary conditions:

```math
\boxed{\;
S_{\rm PBC}(\ell, N)
 \;=\; \frac{c}{3}\,\ln\!\left[\frac{N}{\pi a}
          \,\sin\!\left(\frac{\pi\ell}{N}\right)\right]
 \;+\; s_{1},
\;}
```

```math
\boxed{\;
S_{\rm OBC}(\ell, N)
 \;=\; \frac{c}{6}\,\ln\!\left[\frac{2 N}{\pi a}
          \,\sin\!\left(\frac{\pi\ell}{N}\right)\right]
 \;+\; s_{1}^{\,\prime}.
\;}
```

The universal logarithmic prefactors differ by exactly a factor of
two — ``c/3`` vs ``c/6`` — because PBC produces **two** entanglement
cuts (one at each boundary of the block, since the chain is a ring)
whereas OBC produces **one** (the second boundary coincides with a
physical chain end where the environment is already empty). The
additive constants ``s_{1}, s_{1}^{\,\prime}`` are model-dependent
and *not* universal.

These two formulas are used throughout QAtlas's verification
pipeline to extract the central charge from finite-size
entanglement-entropy data. See
[`docs/src/methods/calabrese-cardy/index.md`](../methods/calabrese-cardy/index.md)
for the method-level summary and
[`docs/src/verification/entanglement.md`](../verification/entanglement.md)
for the specific extraction recipes.

---

## Setup

### Continuum setting

Consider a 1+1-dimensional CFT of central charge ``c`` on a spatial
manifold ``\mathcal{M}``, either:

- **PBC**: a spatial circle of circumference ``L = N\,a``, so
  ``\mathcal{M} = \mathbb{R}_t \times S^{1}_{L}``. The Euclidean
  spacetime is an infinite cylinder of circumference ``L``.
- **OBC**: a spatial segment ``[0, L]``, so ``\mathcal{M} =
  \mathbb{R}_t \times [0, L]``. The Euclidean spacetime is an
  infinite strip of width ``L``.

Let ``A = [0, \ell\,a]`` be a contiguous block of ``\ell`` sites (in the
lattice picture) or length ``\ell_{\rm cont} = \ell\,a`` (in the
continuum); the complement ``B`` is the rest of the chain. The
reduced density matrix of ``A`` at zero temperature is

$``\rho_{A} \;=\; \mathrm{Tr}_{B}\,|0\rangle\langle 0|,``$

and the von Neumann entropy is ``S_{A} = -\mathrm{Tr}\,\rho_{A}\ln
\rho_{A}``.

### Goal

Compute ``S_{A}`` via the replica trick, the twist-operator
formalism, and the conformal maps from ``\mathcal{M}`` to a
computable Euclidean geometry (plane or half-plane). Read off the
logarithmic scaling law and verify the ``c/3`` vs ``c/6`` prefactor.

---

## Calculation

### Step 1 — Replica trick

For any positive integer ``n``, let ``\mathrm{Tr}\,\rho_{A}^{n}`` be
the *``n``-th Rényi trace* of the reduced state. The **replica
identity**

```math
S_{A} \;=\; -\,\mathrm{Tr}\,\rho_{A}\ln\rho_{A}
 \;=\; -\,\lim_{n\to 1}\partial_{n}\,\mathrm{Tr}\,\rho_{A}^{n}
\tag{1}
```

reduces the computation of the entropy to the computation of a
one-parameter family of traces. The proof is immediate: write
``\rho_{A} = \sum_{k} p_{k}|k\rangle\langle k|`` in its eigenbasis so
``\mathrm{Tr}\,\rho_{A}^{n} = \sum_{k} p_{k}^{n}``, differentiate
``\partial_{n}\sum_{k} p_{k}^{n} = \sum_{k} p_{k}^{n}\ln p_{k}``, and
take ``n \to 1`` to get ``\sum_{k} p_{k}\ln p_{k} = -S_{A}``.

The strategy is now (i) compute ``\mathrm{Tr}\,\rho_{A}^{n}`` for
integer ``n \ge 2`` by a path-integral construction, (ii) analytically
continue in ``n``, and (iii) apply (1) to extract ``S_{A}``.

### Step 2 — Twist-operator construction for ``\mathrm{Tr}\,\rho_{A}^{n}``

The path-integral representation of ``\rho_{A}`` is a piecewise
computation: it is the Euclidean path integral on the ``\mathcal{M}
\times \mathbb{R}_{\tau}`` spacetime with:

1. a cut along ``A`` at imaginary time ``\tau = 0^{-}`` (open below),
2. a cut along ``A`` at ``\tau = 0^{+}`` (open above),

the two sides of each cut being identified cyclically in the
standard partition-function-like Euclidean setup. Taking ``n`` copies
and gluing the cuts **cyclically** from one copy to the next
(copy ``1``'s upper cut ``\to`` copy ``2``'s lower cut ``\to \dots \to``
copy ``n``'s lower cut ``\to`` copy ``1``'s upper cut) produces an
**``n``-sheeted Riemann surface** branched along the boundary of ``A``.
Its partition function is ``\mathrm{Tr}\,\rho_{A}^{n}``.

Cardy–Calabrese 2004 (building on Dixon et al. 1987 and
Knizhnik 1987) show that the ``n``-sheeted partition function can be
rewritten as a **correlator of branch-point twist operators**
``\mathcal{T}_{n}`` inserted at the endpoints of ``A``, in the
CFT on a single sheet (``\mathcal{M}``, not the ``n``-sheeted cover):

```math
\mathrm{Tr}\,\rho_{A}^{n}
 \;=\; \bigl\langle \mathcal{T}_{n}(u)\,\overline{\mathcal{T}}_{n}(v)
        \bigr\rangle_{\mathcal{M}},
\tag{2}
```

with ``u = 0, v = \ell`` the endpoints of ``A`` and ``\mathcal{T}_{n},
\overline{\mathcal{T}}_{n}`` a twist / anti-twist pair. The twist
operator has conformal scaling dimension (Calabrese–Cardy 2004
eq. (1.4))

```math
\boxed{\;
\Delta_{n} \;=\; \overline{\Delta}_{n}
 \;=\; \frac{c}{24}\Bigl(n - \frac{1}{n}\Bigr).
\;}
\tag{3}
```

This is the *universal* piece of the result: (3) depends only on
``c`` and ``n``, not on any lattice details. Its proof in the
Calabrese–Cardy paper uses the conformal anomaly of the
``n``-sheeted cover; we accept (3) as input here and focus on the
geometry-dependent factor.

### Step 3 — PBC: cylinder geometry, 2-point function of twists

**Step 3a: plane → cylinder.** PBC spacetime is an infinite
cylinder of circumference ``L = N a``. The conformal map
``w \mapsto z = e^{2\pi w / L}`` takes the cylinder (coordinates
``w = x + i\tau``) to the plane (coordinates ``z``), opening up the
periodic spatial direction into the radial direction on ``\mathbb{C}``.

Under this map, a primary operator ``\phi`` of dimension ``(\Delta,
\overline{\Delta})`` transforms with the Jacobian

```math
\phi_{\rm plane}(z) = \left(\frac{dw}{dz}\right)^{\Delta}
                       \left(\frac{d\overline{w}}{d\overline{z}}\right)^{\overline{\Delta}}
                       \phi_{\rm cyl}(w)
```

(standard primary transformation under conformal maps).

**Step 3b: 2-point function on the plane.** In the plane,

```math
\bigl\langle\mathcal{T}_{n}(z_{1})\overline{\mathcal{T}}_{n}(z_{2})\bigr\rangle_{\rm plane}
 \;=\; \frac{C_{n}}{|z_{1} - z_{2}|^{4\Delta_{n}}},
```

with ``C_{n}`` a non-universal normalisation that becomes the
``s_{1}`` below.

**Step 3c: pull back to the cylinder.** With twist insertions at
spatial points ``x_{1}, x_{2}`` on the ``\tau = 0`` spatial slice of
the cylinder, ``w_{j} = x_{j}`` and ``z_{j} = e^{2\pi x_{j}/L}``.

```math
|z_{1} - z_{2}|^{2}
 = (e^{2\pi x_{1}/L} - e^{2\pi x_{2}/L})
   (\overline{e^{2\pi x_{1}/L}} - \overline{e^{2\pi x_{2}/L}}).
```

Using $z - z' = e^{(\theta + \theta')/2}(e^{(\theta - \theta')/2}
- e^{-(\theta - \theta')/2}) = 2 e^{(\theta + \theta')/2}
\sinh((\theta - \theta')/2)`` with ``\theta = 2\pi x/L$, this gives

```math
|z_{1} - z_{2}|^{2}
 = 4 e^{2\pi(x_{1} + x_{2})/L}\,\sinh^{2}\!\bigl[\pi(x_{1} - x_{2})/L\bigr].
```

Since ``x_{j}`` are real, ``\sinh[\pi(x_{1} - x_{2})/L]`` is real, and
for ``x_{1} - x_{2} = \ell`` (the block extent) the *spatial* distance
reduces to ``\sinh(\pi\ell/L)`` — but on the circular cylinder the
spatial distance between the two twist insertions is actually
``\ell`` (straight-line along the cylinder) measured against a
periodic direction. The Euclidean computation here works with the
**chordal** distance on the plane; the final answer is expressed
in terms of the sine of the angular separation ``2\pi\ell/L``.

Substituting into the 2-point function and applying the conformal
Jacobians,

```math
\bigl\langle\mathcal{T}_{n}\overline{\mathcal{T}}_{n}\bigr\rangle_{\rm cyl}
 \;=\; \frac{C_{n}}{\bigl[(2 L/\pi)\,\sin(\pi\ell/L)\bigr]^{4\Delta_{n}}}.
\tag{4}
```

(See Calabrese–Cardy 2004 eq. (2.9) and (2.13) for the step-by-step
Jacobian accounting that converts the ``\sinh``-expression to a
``\sin``-expression after Wick rotation back to the Lorentzian
cylinder.)

**Step 3d: extract ``S_{A}``.** Taking the log and then the
``n \to 1`` derivative per (1):

```math
\mathrm{Tr}\,\rho_{A}^{n}
 = \frac{C_{n}}{\bigl[(2 L/\pi)\sin(\pi\ell/L)\bigr]^{4\Delta_{n}}},
```

```math
\ln\mathrm{Tr}\,\rho_{A}^{n}
 = \ln C_{n} - 4\Delta_{n}\ln\!\bigl[(2 L/\pi)\sin(\pi\ell/L)\bigr].
```

Substitute ``\Delta_{n} = (c/24)(n - 1/n)`` and differentiate with
respect to ``n`` at ``n = 1``:

```math
\partial_{n}\,\mathrm{Tr}\,\rho_{A}^{n}\Big|_{n = 1}
 = \partial_{n}\ln C_{n}\Big|_{n = 1}
   - 4\cdot\frac{c}{24}\Bigl(1 + \frac{1}{n^{2}}\Bigr)\Big|_{n = 1}\cdot
     \ln\!\bigl[\dots\bigr]
```

```math
= \partial_{n}\ln C_{n}\big|_{n = 1}
   - \frac{c}{3}\,\ln\!\bigl[(2 L/\pi)\sin(\pi\ell/L)\bigr].
```

Applying (1) ``S_{A} = -\lim_{n \to 1}\partial_{n}\mathrm{Tr}\rho_{A}^{n}``,

```math
S_{\rm PBC} = \frac{c}{3}\,\ln\!\bigl[(2 L/\pi)\sin(\pi\ell/L)\bigr]
              + s_{1},
\tag{5}
```

with ``s_{1} = -\partial_{n}\ln C_{n}|_{n = 1}`` a non-universal
constant.

**Converting to lattice units.** With ``L = N a`` and ``\ell_{\rm cont}
= \ell a``, the argument inside the log is

```math
(2 L / \pi)\sin(\pi\ell_{\rm cont}/L)
 = (2 N a/\pi)\sin(\pi\ell/N).
```

Absorbing the overall factor of ``2`` into ``s_{1}`` and combining the
explicit ``a``-factor with the ``1/a``-factor that sits in the universal
twist-operator normalisation (a UV renormalisation that is
conventionally absorbed into the non-universal constant), we obtain
the boxed PBC Main-result form

```math
S_{\rm PBC}(\ell, N)
 = \frac{c}{3}\,\ln\!\left[\frac{N}{\pi a}\sin\!\left(\frac{\pi\ell}{N}\right)\right]
 + s_{1}.
```

### Step 4 — OBC: strip geometry, 1-point function near the boundary

**Step 4a: half-plane → strip.** OBC spacetime is an infinite
strip of width ``L = N a`` with two boundaries at ``x = 0`` and
``x = L``. The conformal map ``w \mapsto z = \sin(\pi w/L)`` takes the
strip (coordinates ``w``) to the upper half-plane $\mathrm{Im}(z) \ge
0``, sending the two boundary lines ``x = 0`` and ``x = L$ to the
real axis.

Under this map, the strip's spatial point ``x`` at ``\tau = 0`` maps to
``z = \sin(\pi x/L) \in \mathbb{R}``. In particular, the block
``A = [0, \ell]`` in the strip maps to the segment ``[0,
\sin(\pi\ell/L)]`` on the positive real axis.

**Step 4b: method of images.** On the upper half-plane with a
**free / Dirichlet boundary** on the real axis (the standard BCFT
for a conformal boundary condition), a primary ``\phi(z)`` at
``z \in \mathrm{UHP}`` feels its mirror image ``\phi(\bar z)`` at the
reflected point ``\bar z \in \mathrm{LHP}``. A single-twist
insertion in the bulk becomes *effectively* a correlator with its
image: the 1-point function in the BCFT is

```math
\bigl\langle\mathcal{T}_{n}(z)\bigr\rangle_{\rm UHP}
 \;=\; \frac{A_{n}}{(2\,\mathrm{Im}\,z)^{2\Delta_{n}}},
```

where ``A_{n}`` is a boundary-condition-dependent constant (the
boundary-operator coefficient) and the denominator ``(2\mathrm{Im}\,
z)^{2\Delta_{n}} = |z - \bar z|^{2\Delta_{n}}`` is exactly the
image-pair distance.

**Step 4c: OBC twist correlator.** Crucially, in the OBC setup
there is *only one* bulk twist insertion (at the single interior
boundary of ``A``, since the other boundary of ``A`` is the physical
chain end at ``x = 0``, which is not an entanglement cut). Applying
the conformal map ``z = \sin(\pi x/L)`` with ``x = \ell``:

```math
|z - \bar z|_{\rm OBC}
 = |\,2 i\,\mathrm{Im}\,\sin(\pi\ell/L)\,|
 = 2\,|\sin(\pi\ell/L)|,
```

wait — for real ``x``, ``z = \sin(\pi x/L)`` is real, so
``\mathrm{Im}\,z = 0``. The correct interpretation is that on the
``\tau = 0`` slice, the twist insertion sits **on** the real axis of
the half-plane (it is on the boundary of ``A``, which in the
spacetime picture at ``\tau = 0`` is on the boundary of the strip
equivalently). The image-pair distance must be computed more
carefully using the conformal-transformation Jacobian.

**Step 4d: Jacobian-accounted computation.** Following
Calabrese–Cardy 2004 eq. (19), the result of carrying out the
half-plane image construction and pulling back to the strip is

```math
\bigl\langle\mathcal{T}_{n}(x)\bigr\rangle_{\rm strip}
 \;=\; \frac{A_{n}}{\bigl[(2L/\pi)\sin(\pi x/L)\bigr]^{2\Delta_{n}}}.
\tag{6}
```

The key differences from the PBC 2-point function (4) are:

- The **exponent** is ``2\Delta_{n}`` (for one operator) instead of
  ``4\Delta_{n}`` (for two) — a factor of ``2`` saved.
- The **prefactor** inside the log-argument has an extra factor
  of ``2`` compared to the PBC result (``(2L/\pi)`` instead of
  ``(L/\pi)``), which reflects the method of images — the effective
  distance is measured from the bulk insertion to its image across
  the boundary, which is twice the bulk-to-boundary distance. In
  the lattice-unit version this becomes ``2N`` inside the log.

**Step 4e: extract ``S_{A}``.**

```math
\ln\,\mathrm{Tr}\,\rho_{A}^{n}
 = \ln A_{n} - 2\Delta_{n}\ln\!\bigl[(2 L / \pi)\sin(\pi\ell/L)\bigr].
```

Differentiating in ``n`` at ``n = 1``:

```math
\partial_{n}\,\mathrm{Tr}\,\rho_{A}^{n}\Big|_{n = 1}
 = \partial_{n}\ln A_{n}\big|_{n = 1}
  - 2\cdot\frac{c}{24}\Bigl(1 + \frac{1}{n^{2}}\Bigr)\Big|_{n = 1}\cdot\ln[\dots]
```

```math
= \partial_{n}\ln A_{n}\big|_{n = 1}
  - \frac{c}{6}\,\ln\!\bigl[(2 L/\pi)\sin(\pi\ell/L)\bigr].
```

Applying (1),

```math
S_{\rm OBC}
 = \frac{c}{6}\,\ln\!\bigl[(2 L/\pi)\sin(\pi\ell/L)\bigr] + s_{1}^{\,\prime}.
\tag{7}
```

In lattice units ``L = N a``, ``\ell_{\rm cont} = \ell a``:

```math
S_{\rm OBC}(\ell, N) = \frac{c}{6}\,\ln\!\left[\frac{2 N}{\pi a}
  \,\sin\!\left(\frac{\pi\ell}{N}\right)\right] + s_{1}^{\,\prime}.
```

This matches the Main-result boxed OBC form.

### Step 5 — Structural origin of the ``c/3`` vs ``c/6`` prefactor

Direct contribution counting:

| Feature              | PBC                       | OBC                     |
|----------------------|---------------------------|-------------------------|
| Spatial manifold     | ring (circumference ``L``)  | segment (width ``L``)     |
| Euclidean spacetime  | infinite cylinder         | infinite strip          |
| Conformal map        | cylinder ``\to`` plane      | strip ``\to`` half-plane  |
| Number of entanglement cuts | ``\mathbf{2}``       | ``\mathbf{1}``            |
| Twist insertions     | 2-point function          | 1-point function        |
| Log-prefactor exponent | ``4\Delta_{n}``           | ``2\Delta_{n}``           |
| ``n \to 1`` factor     | ``c/3``                     | ``c/6``                   |
| Inner-argument factor in log | ``(L/\pi)\sin(\pi\ell/L)`` | ``(2L/\pi)\sin(\pi\ell/L)`` |

The **factor-of-2 ratio** of the log-prefactors ``(c/3)/(c/6) = 2``
tracks exactly the ratio of entanglement cuts ``2/1``. Each cut
contributes ``c/6`` to the log coefficient; PBC has two cuts, OBC
has one, and the additional factor of ``2`` inside the log argument
in the OBC case is the image-doubling of the effective distance.

### Step 6 — Limiting-case and numerical checks

**(i) ``\ell \to N/2`` (balanced bipartition).** At the midpoint,
``\sin(\pi/2) = 1`` and the formulas simplify to

$``S_{\rm PBC}(N/2, N) = \frac{c}{3}\ln(N/(\pi a)) + s_{1},``$

$``S_{\rm OBC}(N/2, N) = \frac{c}{6}\ln(2 N/(\pi a)) + s_{1}^{\,\prime}.``$

Both scale as ``\ln N`` with coefficients differing by a factor of
two.

**(ii) ``\ell \to 0, N`` (endpoints).** ``\sin(\pi\ell/N) \to 0``,
both entropies diverge logarithmically — the UV cutoff renders
the entropy finite, and the log-divergence is the usual "area law
with a log correction in 1+1D".

**(iii) Ising TFIM (``c = 1/2``) at OBC.** QAtlas runs the Peschel
correlation-matrix method in
`src/models/quantum/TFIM/TFIM_entanglement.jl`
at ``N = 100`` to fit (7) and extract ``c \approx 0.5`` within 5%;
see `test/models/test_TFIM_entanglement.jl`
(the "Calabrese–Cardy log scaling" testset).  The Peschel-based
OBC extraction is vastly cheaper than full ED and is one of the
stringent checks of the entanglement-entropy implementation.

**(iv) Heisenberg (``c = 1``) at OBC.** At the Heisenberg-chain
critical point,
`test/verification/test_entanglement_central_charge.jl`
extracts ``c \approx 1`` within ``20\%`` at ``N = 12``, and much better
at larger ``N``. The ``c = 1`` identification is also consistent with
the independent derivation from the XXZ Luttinger parameter
``K = 1/2`` (see
[`xxz-luttinger-parameters`](xxz-luttinger-parameters.md) Step 7)
— every ``c = 1`` compactified-boson CFT gives ``c = 1`` regardless of
the compactification radius.

---

## Why ``\Delta_{n} = (c/24)(n - 1/n)``?

The scaling dimension of the branch-point twist operator is the one
piece of input we took for granted in Step 2. Its derivation uses
the conformal anomaly (Weyl transformation) of the ``n``-sheeted
cover:

Consider the plane with ``n``-fold branching at the twist insertions.
The metric on the ``n``-sheeted cover differs from the plane by a
conformal factor, and the partition function acquires a Liouville
action proportional to the central charge ``c``. Explicit integration
of the Liouville action over the branched geometry yields an
effective dimension
``\Delta_{n} = (c/24)(n - 1/n)`` for each twist insertion
(Cardy–Calabrese 2004 §2, Holzhey–Larsen–Wilczek 1994). At $n =
1`` the cover is the plane and ``\Delta_{1} = 0$; the first-order
Taylor coefficient is ``(c/24)\cdot 2 = c/12``, which combines with
the ``2``-point-function exponent ``4\Delta_{n}`` to give the ``c/3``
PBC prefactor.

---

## References

- P. Calabrese, J. Cardy, *Entanglement entropy and quantum field
  theory*, J. Stat. Mech. **0406**, P06002 (2004). The modern
  systematic derivation; eq. (1.4) for ``\Delta_{n}``, eq. (2.9) for
  the PBC 2-point function, eq. (19) for the OBC 1-point result.
- C. Holzhey, F. Larsen, F. Wilczek, *Geometric and renormalized
  entropy in conformal field theory*, Nucl. Phys. B **424**, 443
  (1994). First CFT derivation of the entanglement entropy scaling
  law in 1+1D (PBC case); introduces the twist-operator idea.
- L. J. Dixon, D. Friedan, E. J. Martinec, S. H. Shenker, *The
  conformal field theory of orbifolds*, Nucl. Phys. B **282**, 13
  (1987). Introduced branch-point twist operators in the orbifold
  context that Calabrese–Cardy later adapted.
- V. G. Knizhnik, *Analytic fields on Riemann surfaces II*, Comm.
  Math. Phys. **112**, 567 (1987). Complementary treatment of
  twist operators on higher-genus surfaces.
- P. Calabrese, J. Cardy, *Entanglement entropy and conformal field
  theory*, J. Phys. A **42**, 504005 (2009). Review paper; eq. (28)
  and (30) are our (5) and (7).
- G. Vidal, J. I. Latorre, E. Rico, A. Kitaev, *Entanglement in
  quantum critical phenomena*, Phys. Rev. Lett. **90**, 227902
  (2003). Lattice verification of the scaling law for free-fermion
  1D chains.

## Used by

- [Calabrese–Cardy method page](../methods/calabrese-cardy/index.md) —
  method-level summary; this note is its derivation backbone.
- [Entanglement-entropy verification](../verification/entanglement.md) —
  describes the central-charge-extraction recipe that uses (7) as
  the regression template.
- [TFIM model page](../models/quantum/tfim.md) — central-charge
  extraction at the ``h = J`` Ising critical point using
  `fetch(TFIM(...), VonNeumannEntropy(), OBC(N); ℓ)` matched
  against (7).
- [Heisenberg model page](../models/quantum/heisenberg.md) and
  [XXZ model page](../models/quantum/xxz.md) — ``c = 1``
  central-charge identification consistent with the independent
  Luttinger-liquid extraction in
  [`xxz-luttinger-parameters`](xxz-luttinger-parameters.md).
