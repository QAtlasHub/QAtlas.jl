# Ising CFT: Primary Operators and Scaling Dimensions

## Main result

The 2D classical Ising model at the self-dual critical point
``\sinh(2\beta_c J) = 1`` is described at long distances by the
simplest non-trivial unitary conformal field theory — the
**Virasoro minimal model ``\mathcal{M}(3, 4)``** of central charge
``c = 1/2``. This CFT contains exactly **three** primary operators,
whose left-moving conformal weights ``h`` (equal to the right-moving
``\bar h`` since the theory is diagonal) and full scaling dimensions
``\Delta = h + \bar h = 2h`` are

```math
\boxed{\;
\begin{array}{c|c|c|c}
\text{Operator} & h = \bar h & \Delta & \text{Lattice correspondence}\\ \hline
\mathbb{I} \text{ (identity)} & 0 & 0 & \text{1} \\
\sigma \text{ (spin)} & 1/16 & 1/8 & \lim_{a\to 0}\, \sigma^{z}(a\,r) \\
\varepsilon \text{ (energy)} & 1/2 & 1 & \lim_{a\to 0}\,\sigma^{z}_{i}\sigma^{z}_{i+1} - \langle \sigma^{z}\sigma^{z}\rangle_{c}\\
\end{array}
\;}
```

The three scaling dimensions ``(0, 1/8, 1)`` are the microscopic
origin of *every* 2D Ising critical exponent (via standard CFT
scaling relations):

```math
\boxed{\;\eta = 2\Delta_{\sigma} = \tfrac{1}{4},\qquad
         \nu = \tfrac{1}{d - \Delta_{\varepsilon}} = 1,\qquad
         \beta = \nu\,\Delta_{\sigma} = \tfrac{1}{8},\qquad
         \alpha = 2 - \nu\,d = 0\;(\log),\qquad
         \delta = \tfrac{d + 2 - \eta}{d - 2 + \eta} = 15.\;}
```

In particular ``\beta = 1/8`` is exactly the Yang-magnetisation
exponent derived independently from the Toeplitz determinant in
[`yang-magnetization-toeplitz`](yang-magnetization-toeplitz.md), and
the ``\eta = 1/4`` correlation-function exponent matches
``2\Delta_{\sigma}`` via Fisher's scaling relation.

The ``(\sigma, \varepsilon)`` pair of non-trivial primaries is
responsible for the two distinct massive integrable perturbations
of the Ising CFT:

- ``\varepsilon`` perturbation → free massive Majorana fermion
  (thermal perturbation, ``T \ne T_{c}``);
- ``\sigma`` perturbation → interacting integrable QFT with ``E_{8}``
  mass spectrum (magnetic perturbation, see
  [`e8-mass-spectrum-derivation`](e8-mass-spectrum-derivation.md)
  and [`ising-cft-magnetic-perturbation`](ising-cft-magnetic-perturbation.md)).

---

## Setup

### Virasoro algebra

A 2D CFT with a unique vacuum is characterised by two copies of
the Virasoro algebra (holomorphic + antiholomorphic),

```math
[L_{m}, L_{n}]
 \;=\; (m - n)\,L_{m + n} \;+\; \frac{c}{12}\,m(m^{2} - 1)\,\delta_{m + n, 0},
```

with the same relation for ``\bar L_{n}`` and $[L_{m}, \bar L_{n}]
= 0``. The central charge ``c$ parametrises the anomaly in the
trace of the stress tensor under conformal rescalings; for the 2D
Ising model ``c = 1/2`` as we confirm in Step 4.

**Primary operators** are operators ``\phi(z, \bar z)`` on which
``L_{n}`` and ``\bar L_{n}`` for ``n > 0`` annihilate the state
``|\phi\rangle = \phi(0)|0\rangle``:

```math
L_{n}|\phi\rangle = 0,\qquad \bar L_{n}|\phi\rangle = 0\qquad
 (n > 0),
```

and ``L_{0}|\phi\rangle = h\,|\phi\rangle``, $\bar L_{0}|\phi\rangle =
\bar h\,|\phi\rangle``. The eigenvalues ``(h, \bar h)$ are the
**conformal weights**; the scaling dimension is ``\Delta = h + \bar h``
and the spin is ``s = h - \bar h``. All scaling operators in a CFT
are either primaries or *descendants* (acted on by ``L_{-n}`` with
``n > 0``).

### Minimal models

BPZ 1984 proved that unitary CFTs with ``c < 1`` are **discrete**:
the only allowed values of the central charge are

```math
c_{p, p'} \;=\; 1 \;-\; \frac{6\,(p' - p)^{2}}{p\,p'},
\qquad p, p' \in \mathbb{Z}_{\ge 2},\ \gcd(p, p') = 1,
\tag{1}
```

with the pair ``(p, p')`` labelling the model
``\mathcal{M}(p, p')``. Each ``\mathcal{M}(p, p')`` has **finitely
many** primary operators, listed by the **Kac table** (Step 3).

### Goal

Derive:
1. ``\mathcal{M}(3, 4)`` has central charge ``c = 1/2``.
2. The Kac table of ``\mathcal{M}(3, 4)`` has exactly three
   inequivalent primary operators.
3. Their conformal weights are ``h \in \{0, 1/16, 1/2\}``.
4. They identify with ``\{\mathbb{I}, \sigma, \varepsilon\}`` of the
   lattice Ising model, reproducing the known critical exponents.

---

## Calculation

### Step 1 — Kac formula for conformal weights

For any minimal model ``\mathcal{M}(p, p')``, the conformal weights
of primary operators are given by the Kac formula

```math
\boxed{\;
h_{r, s} \;=\; \frac{(p'\,r - p\,s)^{2} - (p' - p)^{2}}{4\,p\,p'},
\qquad 1 \le r \le p - 1,\ \ 1 \le s \le p' - 1.
\;}
\tag{2}
```

This is the statement that a primary labelled by ``(r, s)`` has a
null descendant at level ``r\,s`` in its Verma module, which
constrains its four-point correlator to satisfy a linear ODE and
fixes ``h_{r, s}`` uniquely. The derivation requires the Kac
determinant formula (Kac 1979) for the Gram matrix of a Verma
module; see Di Francesco–Mathieu–Sénéchal 1997 §7.3 for the full
treatment.

The Kac table has the **reflection symmetry**

```math
h_{r, s} \;=\; h_{p - r,\, p' - s},
```

which follows from (2) by substituting ``r \to p - r``, $s \to p'
- s$:

```math
(p'(p - r) - p(p' - s))^{2}
 = (p'\,p - p'\,r - p\,p' + p\,s)^{2}
 = (p\,s - p'\,r)^{2}
 = (p'\,r - p\,s)^{2},
```

so the numerator of (2) is invariant and ``h_{r,s} = h_{p-r, p'-s}``.
This symmetry is the reason the Kac table has roughly half as many
distinct entries as the ``(p-1)(p'-1)`` rectangle.

### Step 2 — ``\mathcal{M}(3, 4)``: the Ising CFT

Substitute ``p = 3``, ``p' = 4`` into (1):

```math
c_{3, 4} \;=\; 1 - \frac{6(4 - 3)^{2}}{3\cdot 4}
 \;=\; 1 - \frac{6}{12} \;=\; \frac{1}{2}.
\tag{3}
```

The central charge ``c = 1/2`` is the smallest positive value in the
minimal-model series (the ``c_{2, 3}`` model is trivial with only
the identity operator; the next nontrivial case is $\mathcal{M}(3,
4)$). It is also the central charge of a single free Majorana
fermion — not a coincidence, as we discuss at the end.

The Kac rectangle for ``\mathcal{M}(3, 4)`` has $(p - 1)(p' - 1) =
2 \cdot 3 = 6`` entries labelled by ``(r, s)`` with ``r \in \{1, 2\}$
and ``s \in \{1, 2, 3\}``.

### Step 3 — Kac table for ``\mathcal{M}(3, 4)``

Compute ``h_{r, s}`` for all six ``(r, s)`` pairs using (2) with
``p = 3``, ``p' = 4``, ``p' - p = 1``, ``4 p p' = 48``:

```math
h_{r, s} \;=\; \frac{(4 r - 3 s)^{2} - 1}{48}.
```

Enumerate:

| ``(r, s)`` | ``4r - 3s`` | ``(4r - 3s)^{2} - 1`` | ``h_{r, s}`` |
|----------|-----------|---------------------|------------|
| ``(1, 1)`` | ``1`` | ``0`` | ``0`` |
| ``(1, 2)`` | ``-2`` | ``3`` | ``1/16`` |
| ``(1, 3)`` | ``-5`` | ``24`` | ``1/2`` |
| ``(2, 1)`` | ``5`` | ``24`` | ``1/2`` |
| ``(2, 2)`` | ``2`` | ``3`` | ``1/16`` |
| ``(2, 3)`` | ``-1`` | ``0`` | ``0`` |

Apply the Kac-symmetry ``h_{r, s} = h_{p - r, p' - s}``:

- ``(1, 1) \leftrightarrow (2, 3)``: both ``h = 0``. ✓
- ``(1, 2) \leftrightarrow (2, 2)``: both ``h = 1/16``. ✓
- ``(1, 3) \leftrightarrow (2, 1)``: both ``h = 1/2``. ✓

So the six entries collapse to **three inequivalent conformal
weights**

```math
\boxed{\;h \in \{0,\; 1/16,\; 1/2\}.\;}
\tag{4}
```

Each of these is a self-conjugate primary (for a diagonal CFT like
the Ising model, ``(h, \bar h)`` are equal — left = right).
Therefore the 2D Ising CFT has exactly **three scalar primary
operators**, with full scaling dimensions

```math
\Delta \in \{0,\ 1/8,\ 1\}.
\tag{5}
```

### Step 4 — Identification with lattice operators

We now identify the three primaries with physical operators on the
Ising lattice. The identification uses the **scaling dimension** of
lattice operators in the critical two-point function:

```math
\bigl\langle\mathcal{O}(\mathbf{r})\,\mathcal{O}(\mathbf{0})\bigr\rangle_{T_{c}}
 \;\sim\; |\mathbf{r}|^{-2\Delta_{\mathcal{O}}}.
\tag{6}
```

**Identity ``\mathbb{I}``, ``\Delta = 0``**: the trivial constant
operator. Its correlator is the vacuum expectation, constant in
``|\mathbf{r}|``.

**Spin ``\sigma``, ``\Delta = 1/8``**: the continuum limit of the Ising
spin variable ``\sigma^{z}_{i}``. The two-point function decays as

```math
\langle\sigma^{z}(\mathbf{r})\,\sigma^{z}(\mathbf{0})\rangle_{T_{c}}
 \;\sim\; |\mathbf{r}|^{-2\Delta_{\sigma}}
 \;=\; |\mathbf{r}|^{-1/4}.
```

This matches the ``\eta = 1/4`` Fisher exponent of the 2D Ising
model (Fisher 1964; Stanley 1971 §12.4), derivable from the
Toeplitz-determinant analysis of
[`yang-magnetization-toeplitz`](yang-magnetization-toeplitz.md).

**Energy ``\varepsilon``, ``\Delta = 1``**: the continuum limit of
the bond-energy operator
$\varepsilon_{i} \propto \sigma^{z}_{i}\sigma^{z}_{i+1} -
\langle\sigma^{z}\sigma^{z}\rangle_{c}$ (the connected bond-bond
correlator). Its conformal dimension ``1`` reflects the thermal
operator — perturbing by ``\varepsilon`` tunes ``T`` away from ``T_{c}``,
with ``T - T_{c} \propto \int d^{2}x\,\varepsilon(x)``. The relation
``\nu = 1/(d - \Delta_{\varepsilon}) = 1/(2 - 1) = 1`` reads off the
correlation-length exponent directly.

The lattice → CFT identifications are

```math
\sigma^{z}_{i} \;\xrightarrow[a \to 0]{}\; a^{\Delta_{\sigma}}\,\sigma(\mathbf{r}_{i})
 \;=\; a^{1/8}\,\sigma(\mathbf{r}_{i}),
```

```math
\sigma^{z}_{i}\sigma^{z}_{i+1} - \langle\sigma^{z}\sigma^{z}\rangle
 \;\xrightarrow[a \to 0]{}\; a^{\Delta_{\varepsilon}}\,\varepsilon(\mathbf{r}_{i})
 \;=\; a^{1}\,\varepsilon(\mathbf{r}_{i}),
```

where the powers of the lattice spacing ``a`` accompany the RG
blocking of each operator under a rescaling ``\mathbf{r} \to b\mathbf{r}``.

### Step 5 — Critical exponents from ``(\Delta_{\sigma}, \Delta_{\varepsilon})``

All critical exponents of the 2D Ising universality class follow
from the two non-trivial scaling dimensions ``\Delta_{\sigma} = 1/8``
and ``\Delta_{\varepsilon} = 1`` via standard CFT / scaling
relations (Cardy 1996 Ch. 3).

**Fisher exponent ``\eta``** (anomalous dimension of the order
parameter). Definition: $\langle\sigma\sigma\rangle \sim
|\mathbf{r}|^{-(d - 2 + \eta)}$. Compare with (6):
``2\Delta_{\sigma} = d - 2 + \eta``, hence

```math
\eta \;=\; 2\Delta_{\sigma} - d + 2
 \;=\; 2\cdot\tfrac{1}{8} - 2 + 2 \;=\; \tfrac{1}{4}.
```

**Correlation-length exponent ``\nu``.** Definition: $\xi \sim
|T - T_{c}|^{-\nu}$. Scaling dimension of the thermal coupling
``\delta K = \beta(T - T_{c})/J`` is ``d - \Delta_{\varepsilon}``
(since ``\int d^{2}x\,\varepsilon`` must be dimensionless). Hence

```math
\nu \;=\; \frac{1}{d - \Delta_{\varepsilon}}
 \;=\; \frac{1}{2 - 1} \;=\; 1.
```

**Order-parameter exponent ``\beta``.** Definition: $M \sim
|T_{c} - T|^{\beta}``. Scaling of ``M`` as ``\xi^{-\Delta_{\sigma}}$
gives ``M \sim \xi^{-\Delta_{\sigma}} \sim |T_{c} - T|^{\nu\Delta_{\sigma}}``,
hence

```math
\beta \;=\; \nu\,\Delta_{\sigma} \;=\; 1\cdot\tfrac{1}{8} \;=\; \tfrac{1}{8}.
```

Cross-check with
[`yang-magnetization-toeplitz`](yang-magnetization-toeplitz.md) Step
7, which obtains ``\beta = 1/8`` from the Szegő-theorem evaluation
of the Toeplitz determinant: both derivations agree to first
principles.

**Susceptibility exponent ``\gamma``.** From Widom scaling $\gamma =
\nu(2 - \eta) = 1 \cdot (2 - 1/4) = 7/4$. Or directly:
the susceptibility scales as ``\int d^{2}x\,\langle\sigma\sigma\rangle``
which diverges as ``\xi^{d - 2\Delta_{\sigma}}``, giving
$\chi \sim |T - T_{c}|^{-\nu(d - 2\Delta_{\sigma})} = |T -
T_{c}|^{-7/4}$.

**Specific-heat exponent ``\alpha``.** Josephson scaling $\alpha = 2
- \nu d = 2 - 2 = 0``, i.e. ``\alpha = 0$ with a logarithmic
divergence (Onsager 1944). This matches the exact 2D Ising specific
heat ``C \sim |\ln|T - T_{c}||``.

**Critical isotherm ``\delta``.** Widom relation $\delta = (d + 2
- \eta)/(d - 2 + \eta) = (4 - 1/4)/(0 + 1/4) = (15/4)/(1/4) = 15$.

Summary of the five independent 2D Ising exponents — all derived
from the Kac table (2) + scaling relations:

| Exponent | CFT formula | Value |
|----------|-------------|-------|
| ``\eta`` | ``2\Delta_{\sigma} - d + 2`` | ``1/4`` |
| ``\nu`` | ``1/(d - \Delta_{\varepsilon})`` | ``1`` |
| ``\beta`` | ``\nu\,\Delta_{\sigma}`` | ``1/8`` |
| ``\gamma`` | ``\nu(2 - \eta) = \nu(d - 2\Delta_{\sigma})`` | ``7/4`` |
| ``\alpha`` | ``2 - \nu d`` | ``0`` (log) |
| ``\delta`` | ``(d + 2 - \eta)/(d - 2 + \eta)`` | ``15`` |

Further scaling-relation cross-checks in
[`ising-scaling-relations`](ising-scaling-relations.md); the QAtlas
stored values are in
[Ising universality class page](../universalities/ising.md).

### Step 6 — Why ``c = 1/2``? The free-Majorana-fermion
construction

The central charge ``c = 1/2`` is also that of a single free Majorana
fermion (``c = N/2`` for ``N`` real Majoranas; a complex Dirac fermion
has ``c = 1``). This is not a coincidence: the 2D Ising CFT **is**
a free Majorana fermion, via the Jordan–Wigner transformation of
the lattice TFIM at the critical point ``h = J``.

From [`jw-tfim-bdg`](jw-tfim-bdg.md) Step 7(iii), the
thermodynamic-limit dispersion of the TFIM at ``h = J`` is
``\Lambda(k) = 4 J\,|\sin(k/2)|``, vanishing linearly at ``k = 0``
(and ``k = 2\pi``). Linearising near ``k = 0`` gives two chiral modes
``\Psi_{R}(x), \Psi_{L}(x)`` forming a single Majorana fermion
doublet. The CFT of the long-wavelength theory is therefore a
*single* free Majorana — central charge ``c = 1/2``. The three
primaries ``\mathbb{I}, \sigma, \varepsilon`` correspond to the
identity, the **disorder operator** (non-local in fermion
language, local in the ``\mathbb{Z}_{2}``-dual Ising formulation),
and the **fermion bilinear** ``\bar\psi\psi`` (= ``\sigma^{z}\sigma^{z}``
after JW) respectively.

The spin–disorder duality between the ``\sigma`` and the dual-theory
spin is the CFT avatar of Kramers–Wannier duality; see
[`kramers-wannier-duality`](kramers-wannier-duality.md).

### Step 7 — Limiting-case and consistency checks

**(i) Identity operator.** ``\Delta_{\mathbb{I}} = 0`` as required of
any CFT: ``\langle\mathbb{I}\rangle = 1`` is constant in ``\mathbf{r}``.

**(ii) Unitarity of ``\mathcal{M}(3, 4)``.** The FQS theorem
(Friedan–Qiu–Shenker 1984) classifies unitary minimal models as
those with ``p' = p + 1``; checking ``\mathcal{M}(3, 4)`` has
``p' - p = 1``, so it is unitary. The three primaries have
non-negative conformal weights (0, 1/16, 1/2), consistent.

**(iii) Numerical verification against finite-size ED.**
``\sigma``-type operator: the two-point function of ``\sigma^{z}`` on
the critical TFIM decays as ``|\mathbf{r}|^{-1/4}``; the exponent is
extracted from QAtlas's closed-form correlator in the
`"criticality: scaling toward CFT exponent -1/4"` testset of
`test/models/test_TFIM_dynamics.jl`, which asserts the effective
doubling-ratio slope at ``N \in \{80, 160, 240\}`` converges
monotonically to ``-1/4`` within ``0.10`` at the largest ``N``.
The central charge ``c = 1/2`` itself is extracted to ``\lesssim 1\%``
from the PBC entanglement entropy in
`test/verification/test_entanglement_central_charge.jl`.

**(iv) Yang-magnetisation exponent ``\beta = 1/8``.** Derived
independently from the Toeplitz determinant in
[`yang-magnetization-toeplitz`](yang-magnetization-toeplitz.md).
Both derivations yield ``\beta = 1/8`` — the microscopic confirmation
of the CFT prediction $\beta = \nu\,\Delta_{\sigma} = 1 \cdot
\tfrac{1}{8}$.

**(v) Logarithmic specific heat.** Onsager 1944 exactly solved the
2D Ising free energy and found ``C \sim -\ln|T - T_{c}|`` near
criticality. The CFT value ``\alpha = 0`` corresponds to this
logarithmic (not power-law) divergence; the ``\alpha = 0`` value is
"degenerate" with a ``\ln`` correction, a common feature of
2D CFTs.

---

## References

- A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov, *Infinite
  conformal symmetry in two-dimensional quantum field theory*,
  Nucl. Phys. B **241**, 333 (1984). Original introduction of
  minimal models; central-charge formula (1) is BPZ eq. (3.24).
- V. G. Kac, *Contravariant form for infinite-dimensional Lie
  algebras and superalgebras*, in *Group Theoretical Methods in
  Physics*, Springer (1979) 441–445. The Kac determinant formula
  that determines the null-vector structure and hence (2).
- D. Friedan, Z. Qiu, S. Shenker, *Conformal invariance,
  unitarity and critical exponents in two dimensions*, Phys. Rev.
  Lett. **52**, 1575 (1984). FQS classification of unitary minimal
  models.
- J. L. Cardy, *Scaling and Renormalization in Statistical
  Physics*, Cambridge University Press (1996), Ch. 3. Critical
  exponents from CFT + scaling relations.
- P. Di Francesco, P. Mathieu, D. Sénéchal, *Conformal Field
  Theory*, Springer (1997), Ch. 7. Pedagogical treatment of
  minimal models; the Kac formula derivation is in §7.3, and the
  Ising-CFT / free-Majorana identification is in §12.
- J. L. Cardy, *Conformal field theory and statistical mechanics*,
  Les Houches lecture notes (2008), arXiv:0807.3472. Modern review.

## Used by

- [Ising universality class](../universalities/ising.md) — the
  stored exponents $\beta = 1/8, \gamma = 7/4, \nu = 1, \eta = 1/4,
  \alpha = 0, \delta = 15$ are read off from the scaling dimensions
  derived here.
- [`yang-magnetization-toeplitz`](yang-magnetization-toeplitz.md) —
  independent microscopic derivation of ``\beta = 1/8`` via the
  Szegő / Toeplitz determinant; agrees with the CFT
  ``\beta = \nu\,\Delta_{\sigma}`` prediction.
- [`ising-scaling-relations`](ising-scaling-relations.md) — the
  scaling relations linking all six exponents are verified in
  Step 5 above.
- [`e8-mass-spectrum-derivation`](e8-mass-spectrum-derivation.md) —
  the ``\sigma``-perturbation channel produces the ``E_{8}`` massive
  integrable theory.
- [`ising-cft-magnetic-perturbation`](ising-cft-magnetic-perturbation.md) —
  conserved-charges mechanism for the ``E_{8}`` perturbation.
- [`kramers-wannier-duality`](kramers-wannier-duality.md) — the
  CFT version of KW duality is the spin↔disorder duality between
  ``\sigma`` and its ``\mathbb{Z}_{2}``-dual primary.
- [TFIM model page](../models/quantum/tfim.md) — ``h = J`` is the
  lattice realisation of the Ising CFT critical point.
