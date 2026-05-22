# Ising Critical Exponents: Scaling Relations

## Main result

The six standard critical exponents $\alpha, \beta, \gamma, \delta,
\nu, \eta$ of a second-order phase transition are not independent.
Given the **Widom scaling hypothesis** that the singular part of the
free-energy density is a generalised homogeneous function of the
reduced temperature and magnetic field,

```math
f_{\rm sing}(t, h) \;=\; |t|^{2 - \alpha}\,\Phi\!\left(\frac{h}{|t|^{\Delta_{h}}}\right),
\tag{0}
```

four exact scaling relations follow:

```math
\boxed{\;
\begin{aligned}
\alpha + 2\beta + \gamma &= 2 & &\text{(Rushbrooke, 1963)}\\
\gamma &= \beta(\delta - 1) & &\text{(Widom, 1965)}\\
\gamma &= \nu(2 - \eta) & &\text{(Fisher, 1969)}\\
2 - \alpha &= d\,\nu & &\text{(Josephson, 1967; valid for }d \le 4\text{)}
\end{aligned}
\;}
```

Four relations among six exponents ``\Rightarrow`` **two independent
exponents**. For the 2D Ising universality class the two
independent exponents can be taken to be the pair
``(\Delta_{\sigma}, \Delta_{\varepsilon}) = (1/8, 1)`` of CFT scaling
dimensions (see
[`ising-cft-primary-operators`](ising-cft-primary-operators.md)),
and all six thermodynamic exponents follow:

```math
\alpha = 0\ (\log),\quad \beta = \tfrac{1}{8},\quad
\gamma = \tfrac{7}{4},\quad \delta = 15,\quad \nu = 1,\quad \eta = \tfrac{1}{4}.
```

All four scaling relations above are satisfied exactly by this
tuple — the QAtlas test suite verifies this using `Rational{Int}`
arithmetic (no floating-point slack) in
`test/standalone/test_universality_exponents.jl`.

---

## Setup

### Six critical exponents

Define six critical exponents in terms of the singular behaviour
of thermodynamic quantities near a second-order phase transition at
``T = T_{c}``, with reduced temperature ``t = (T - T_{c})/T_{c}`` and
applied field ``h``:

| Exponent | Definition | Quantity |
|----------|------------|----------|
| ``\alpha`` | ``C \sim |t|^{-\alpha}`` (``h = 0``) | specific heat (per unit volume) |
| ``\beta`` | ``M \sim (-t)^{\beta}`` (``t < 0,\ h = 0``) | order parameter (spontaneous magnetisation) |
| ``\gamma`` | ``\chi \sim |t|^{-\gamma}`` (``h = 0``) | susceptibility ``\chi = (\partial M/\partial h)_{h = 0}`` |
| ``\delta`` | ``M \sim h^{1/\delta}\,\mathrm{sgn}(h)`` (``t = 0``) | critical isotherm |
| ``\nu`` | ``\xi \sim |t|^{-\nu}`` (``h = 0``) | correlation length |
| ``\eta`` | ``G(\mathbf{r}) \sim |\mathbf{r}|^{-(d - 2 + \eta)}`` (``t = 0``) | anomalous dimension |

The full set ``\{\alpha, \beta, \gamma, \delta, \nu, \eta\}`` looks
*a priori* independent, but the scaling hypothesis forces four
algebraic relations among them.

### Widom scaling hypothesis

**Hypothesis** (Widom 1965). Near ``T_{c}``, the singular part of
the free-energy density has the generalised homogeneous form

```math
f_{\rm sing}(t, h) \;=\; |t|^{2 - \alpha}\,
   \Phi\!\bigl(h\,|t|^{-\Delta_{h}}\bigr),
\tag{1}
```

with ``\Phi(\cdot)`` a smooth scaling function and a single gap
exponent ``\Delta_{h}`` fixed by the two relations we derive below.

The **two-parameter scaling** (1) is the microscopic statement that
``(t, h)`` renormalise as ``(b^{y_{t}} t, b^{y_{h}} h)`` under a
rescaling ``\mathbf{r} \to b\mathbf{r}``, with exponents ``y_{t}, y_{h}``
determined by the CFT scaling dimensions of the thermal (``\varepsilon``)
and magnetic (``\sigma``) relevant operators:

```math
y_{t} = d - \Delta_{\varepsilon},\qquad y_{h} = d - \Delta_{\sigma}.
\tag{2}
```

(The scaling dimension of a relevant coupling is $y = d -
\Delta_{\rm operator}$ — see
[`ising-cft-primary-operators`](ising-cft-primary-operators.md)
Step 5.) In 2D Ising, ``\Delta_{\varepsilon} = 1`` gives ``y_{t} = 1``
and ``\Delta_{\sigma} = 1/8`` gives ``y_{h} = 15/8``.

The scaling form (1) captures the consequences of this
renormalisation-group structure without assuming the specific
microscopic operator content.

### Goal

Derive the four scaling relations (Rushbrooke, Widom, Fisher,
Josephson) from (1) by differentiating the singular free energy and
matching the power-law behaviours of the resulting thermodynamic
quantities.

---

## Calculation

### Step 1 — Consequences of the scaling form (1) for ``M, \chi, C``

From the thermodynamic definitions:

**Order parameter.** ``M = -(\partial f/\partial h)_{t}``. At ``h = 0``,
differentiating (1) once,

```math
M(t, 0) \;=\; -\,|t|^{2 - \alpha - \Delta_{h}}\,\Phi'(0)\cdot
 \mathrm{sgn}(t),
```

vanishing at ``t = 0`` and scaling as ``|t|^{2 - \alpha - \Delta_{h}}``
for ``t < 0`` (below ``T_c``). Comparing to the definition $M \sim
(-t)^{\beta}$:

```math
\boxed{\;\beta \;=\; 2 - \alpha - \Delta_{h}.\;}
\tag{3}
```

**Susceptibility.** ``\chi = -(\partial^{2} f/\partial h^{2})_{t}`` at
``h = 0``:

```math
\chi(t, 0) \;=\; -\,|t|^{2 - \alpha - 2\Delta_{h}}\,\Phi''(0).
```

So

```math
|t|^{-\gamma} \;\sim\; |t|^{2 - \alpha - 2\Delta_{h}}
\quad\Longrightarrow\quad
\boxed{\;\gamma \;=\; 2\Delta_{h} + \alpha - 2.\;}
\tag{4}
```

**Specific heat.** ``C = -T(\partial^{2} f/\partial t^{2})_{h}`` at
``h = 0``:

```math
C(t, 0) \;\sim\; |t|^{-\alpha},
```

matching the definition by construction — the exponent ``\alpha`` in
(1) is precisely the specific-heat exponent (this is why the
``2 - \alpha`` prefactor is used in the scaling form).

**Critical isotherm.** At ``t = 0``, the scaling variable $h|t|^{-\Delta_{h}}
\to \infty``, and ``\Phi(x) \sim x^{p}`` for some power ``p$ as
``x \to \infty``. Requiring ``f`` to be finite at ``t = 0`` fixes
``p = (2 - \alpha)/\Delta_{h}``, so

```math
f(0, h) \;\sim\; h^{(2 - \alpha)/\Delta_{h}}.
```

Then $M(0, h) = -(\partial f/\partial h)_{t = 0} \sim
h^{(2 - \alpha)/\Delta_{h} - 1}``, comparing to ``M \sim h^{1/\delta}$:

```math
\frac{1}{\delta}
 \;=\; \frac{2 - \alpha}{\Delta_{h}} - 1
 \;=\; \frac{2 - \alpha - \Delta_{h}}{\Delta_{h}}
 \;=\; \frac{\beta}{\Delta_{h}},
\quad\Longrightarrow\quad
\boxed{\;\Delta_{h} = \beta\,\delta.\;}
\tag{5}
```

The gap exponent ``\Delta_{h}`` is thus determined by the product
``\beta\,\delta``, sometimes called the "gap" between the order-
parameter and critical-isotherm exponents.

### Step 2 — Rushbrooke: ``\alpha + 2\beta + \gamma = 2``

Add (3) and (4):

```math
\beta + \gamma
 = (2 - \alpha - \Delta_{h}) + (2\Delta_{h} + \alpha - 2)
 = \Delta_{h}.
```

Using (5) ``\Delta_{h} = \beta\,\delta``:

```math
\beta + \gamma \;=\; \beta\,\delta,
\quad\Longrightarrow\quad
\gamma \;=\; \beta(\delta - 1).
```

That is the **Widom relation** (we'll return to it). Meanwhile,
re-arranging (3) directly gives

```math
\Delta_{h} \;=\; 2 - \alpha - \beta.
```

Substitute into (4) $\gamma = 2\Delta_{h} + \alpha - 2 = 2(2 -
\alpha - \beta) + \alpha - 2 = 2 - \alpha - 2\beta$, i.e.

```math
\boxed{\;\alpha + 2\beta + \gamma \;=\; 2.\;}
\tag{RUSHBROOKE}
```

**2D Ising check**: $0 + 2\cdot\tfrac{1}{8} + \tfrac{7}{4}
= 0 + \tfrac{1}{4} + \tfrac{7}{4} = 2$. ✓

Historical note: Rushbrooke 1963 originally proved the *inequality*
``\alpha + 2\beta + \gamma \ge 2`` from thermodynamic stability; the
scaling hypothesis saturates the bound to an equality.

### Step 3 — Widom: ``\gamma = \beta(\delta - 1)``

From the derivation above, the relation $\beta + \gamma =
\beta\,\delta`` is just a rewriting of ``\gamma = \beta(\delta - 1)$:

```math
\boxed{\;\gamma \;=\; \beta(\delta - 1).\;}
\tag{WIDOM}
```

**2D Ising check**: $\tfrac{1}{8}\cdot(15 - 1) = \tfrac{14}{8} =
\tfrac{7}{4}$. ✓

**Alternative derivation from equation of state.** On the critical
isotherm ``t = 0``, the equation of state reads ``M \sim h^{1/\delta}``
from (5). Off the critical isotherm at ``t < 0``, integrating the
``\chi \sim |t|^{-\gamma}`` over ``h \sim 0`` gives
``M \sim \chi\,h \sim |t|^{-\gamma}\,h`` for small ``h``. Matching
the small-``h`` and large-``h`` regimes through the scaling form fixes
``\gamma = \beta(\delta - 1)``; the detailed derivation uses the
short-distance behaviour of ``\Phi``.

### Step 4 — Fisher: ``\gamma = \nu(2 - \eta)``

The susceptibility is the integrated connected correlator

```math
\chi \;=\; \int d^{d}\mathbf{r}\,\bigl\langle s(\mathbf{r})\,
                                            s(\mathbf{0})\bigr\rangle_{c}.
```

Near ``T_{c}`` at ``h = 0``, the correlator has the scaling form

```math
\bigl\langle s(\mathbf{r})\,s(\mathbf{0})\bigr\rangle_{c}
 \;=\; |\mathbf{r}|^{-(d - 2 + \eta)}\,\Psi(|\mathbf{r}|/\xi),
```

with ``\Psi`` a cutoff function that decays rapidly for $|\mathbf{r}|
\gg \xi$. Substituting and switching to the dimensionless variable
``u = |\mathbf{r}|/\xi``,

```math
\chi \;=\; \int\,|\mathbf{r}|^{-(d - 2 + \eta)}\,\Psi(|\mathbf{r}|/\xi)
             \,d^{d}\mathbf{r}
       \;=\; \xi^{2 - \eta}\int u^{-(d - 2 + \eta)}\,\Psi(u)\,
       u^{d - 1}\,du
```
```math
= \xi^{2 - \eta}\int u^{1 - \eta}\,\Psi(u)\,du \;\sim\; \xi^{2 - \eta}.
```

Using ``\xi \sim |t|^{-\nu}``,

```math
\chi \sim |t|^{-\nu(2 - \eta)},
```

and comparing to ``\chi \sim |t|^{-\gamma}``:

```math
\boxed{\;\gamma \;=\; \nu(2 - \eta).\;}
\tag{FISHER}
```

**2D Ising check**: ``1\cdot(2 - 1/4) = 7/4``. ✓

### Step 5 — Josephson (hyperscaling): ``2 - \alpha = d\,\nu``

The singular part of the free-energy density scales as the inverse
of the correlation volume times a characteristic energy:

```math
f_{\rm sing}(t, 0)
 \;\sim\; \frac{k_{B} T}{\xi^{d}(t)}
 \;\sim\; \xi^{-d}
 \;\sim\; |t|^{d\nu}.
```

This is the "hyperscaling" assumption: each correlation volume
contributes order ``k_{B}T`` of free energy, and there are
``\xi^{-d}`` such volumes per unit volume.

On the other hand, from (1) at ``h = 0`` the scaling form gives
``f_{\rm sing}(t, 0) \sim |t|^{2 - \alpha}``. Matching the two:

```math
|t|^{2 - \alpha} \;\sim\; |t|^{d\nu}
\quad\Longrightarrow\quad
\boxed{\;2 - \alpha \;=\; d\,\nu.\;}
\tag{JOSEPHSON}
```

**2D Ising check** (``d = 2``): ``2 - 0 = 2 \cdot 1 = 2``. ✓

**Validity regime.** The hyperscaling argument assumes that
``\xi^{-d}`` is the *only* relevant IR cutoff on the free energy.
This fails above the **upper critical dimension** ``d_{c} = 4`` (for
the Ising universality class), because at ``d > 4`` mean-field
exponents take over and the "hyperscaling" relation $2 - \alpha =
d\nu`` would predict ``2 - 0 = d/2``, i.e. ``d = 4$ — any higher and
the relation breaks. Above ``d_{c}`` the mean-field values
$\alpha = 0, \beta = 1/2, \gamma = 1, \delta = 3, \nu = 1/2,
\eta = 0$ satisfy (RUSHBROOKE), (WIDOM), (FISHER) but violate
(JOSEPHSON) in ``d > 4``.

### Step 6 — All six exponents from ``(\Delta_{\sigma}, \Delta_{\varepsilon})``

From
[`ising-cft-primary-operators`](ising-cft-primary-operators.md) Step
5, the CFT scaling dimensions give:

```math
\eta = 2\Delta_{\sigma} + 2 - d,\qquad
\nu = \frac{1}{d - \Delta_{\varepsilon}},\qquad
\beta = \nu\,\Delta_{\sigma}.
```

Using Widom and Rushbrooke, the remaining three exponents are

```math
\gamma = \beta(\delta - 1),\qquad
\alpha = 2 - 2\beta - \gamma,\qquad
\delta = \frac{d + 2 - \eta}{d - 2 + \eta}.
```

For 2D Ising ($d = 2, \Delta_{\sigma} = 1/8, \Delta_{\varepsilon}
= 1$):

- ``\eta = 2\cdot 1/8 + 2 - 2 = 1/4``
- ``\nu = 1/(2 - 1) = 1``
- ``\beta = 1\cdot 1/8 = 1/8``
- ``\delta = (2 + 2 - 1/4)/(2 - 2 + 1/4) = (15/4)/(1/4) = 15``
- ``\gamma = (1/8)\cdot 14 = 7/4``
- ``\alpha = 2 - 2\cdot 1/8 - 7/4 = 2 - 1/4 - 7/4 = 0`` (log)

All six 2D Ising exponents are determined by the two CFT scaling
dimensions ``(1/8, 1)``.  This is the minimal-data statement of the
2D Ising universality class.

### Step 7 — Algebraic verification (`Rational{Int}`)

The four scaling relations, the six exponents, and their pairwise
consistency can all be verified in exact arithmetic — QAtlas does
this in `test/standalone/test_universality_exponents.jl`
using `Rational{Int}`:

```julia
using Test

(α, β, γ, δ, ν, η) = (0//1, 1//8, 7//4, 15//1, 1//1, 1//4)
d = 2

@test α + 2β + γ == 2            # Rushbrooke
@test γ == β * (δ - 1)           # Widom
@test γ == ν * (2 - η)           # Fisher
@test 2 - α == d * ν             # Josephson
```

All four assertions pass exactly (no floating-point tolerance
needed). The same test is extended in QAtlas to the other
universality classes stored in `src/universalities/*.jl`
(percolation, Potts-``q``, KPZ, ``O(N)`` models) with the appropriate
exponent tuples.

### Step 8 — Limiting-case consistency checks

**(i) Mean-field exponents (``d \ge d_{c} = 4``).** The Landau-theory
values

```math
(\alpha, \beta, \gamma, \delta, \nu, \eta)_{\rm MF} \;=\;
 (0,\, \tfrac{1}{2},\, 1,\, 3,\, \tfrac{1}{2},\, 0)
```

satisfy:
- Rushbrooke: ``0 + 1 + 1 = 2`` ✓
- Widom: ``1 = \tfrac{1}{2}\cdot 2`` ✓
- Fisher: ``1 = \tfrac{1}{2}\cdot 2`` ✓
- Josephson: ``2 - 0 = 4\cdot\tfrac{1}{2} = 2`` ✓ at ``d = 4``; fails
  for ``d > 4`` (predicts ``d\nu = d/2 \ne 2``).

Above ``d_{c} = 4``, Josephson requires ``d\nu = 2`` with ``\nu = 1/2``,
so only ``d = 4``. For ``d > 4`` one must drop the hyperscaling
assumption, replacing ``\xi^{-d}`` with ``\xi^{-d_{c}}`` in the
free-energy argument — a classical *dangerous irrelevant variable*
effect (Fisher 1983).

**(ii) 3D Ising (numerical bootstrap / Monte Carlo).** With
``d = 3`` and the conformal-bootstrap exponents (Kos–Poland–Simmons-
Duffin 2016, Phys. Rev. D **93**, 036404):
``\beta = 0.32642``, ``\gamma = 1.23708``, ``\nu = 0.62997``,
``\eta = 0.03631``, ``\alpha \approx 0.110``, ``\delta \approx 4.79``.

Rushbrooke: ``0.110 + 2\cdot 0.32642 + 1.23708 = 1.99992 \approx 2``. ✓
Widom: ``0.32642\cdot(4.79 - 1) \approx 1.237`` vs ``\gamma = 1.237``. ✓
Fisher: ``0.62997\cdot(2 - 0.03631) = 0.62997\cdot 1.96369 = 1.23708``
vs ``\gamma = 1.23708``. ✓
Josephson: ``2 - 0.110 = 1.890`` vs ``3\cdot 0.62997 = 1.890``. ✓

All four relations hold to the accuracy of the bootstrap
determination (5–6 significant figures).

**(iii) Scaling-relation completeness.** Given the four relations
above, any two exponents determine the other four. The "canonical"
choices are either ``(\alpha, \nu)`` (thermodynamic:
specific-heat + correlation-length), or $(\Delta_{\sigma},
\Delta_{\varepsilon})$ (CFT: magnetic + thermal scaling dimensions).
The latter is more fundamental because $\Delta_{\sigma},
\Delta_{\varepsilon}$ are the defining scaling dimensions of the
RG fixed-point operator algebra, whereas ``\alpha, \nu`` are derived
measurable quantities.

---

## References

- B. Widom, *Equation of state in the neighborhood of the critical
  point*, J. Chem. Phys. **43**, 3898 (1965). Original scaling-form
  hypothesis (1) and derivation of (WIDOM).
- G. S. Rushbrooke, *On the thermodynamics of the critical region
  for the Ising problem*, J. Chem. Phys. **39**, 842 (1963).
  Thermodynamic-stability *inequality* $\alpha + 2\beta + \gamma
  \ge 2$; saturated to equality by scaling.
- M. E. Fisher, *Rigorous inequalities for critical-point
  correlation exponents*, Phys. Rev. **180**, 594 (1969).
  Derivation of (FISHER) from integrated correlator.
- B. D. Josephson, *Inequality for the specific heat. I. Derivation*,
  Proc. Phys. Soc. **92**, 269 (1967) and II. *Application to
  critical phenomena*, ibid. 276. Hyperscaling inequality
  ``2 - \alpha \ge d\nu``; saturated to equality by scaling.
- L. P. Kadanoff, W. Götze, D. Hamblen, R. Hecht, E. A. S. Lewis,
  V. V. Palciauskas, M. Rayl, J. Swift, D. Aspnes, J. Kane,
  *Static phenomena near critical points: theory and experiment*,
  Rev. Mod. Phys. **39**, 395 (1967). Comprehensive review.
- M. E. Fisher, *Scaling, universality and renormalization group
  theory*, in *Critical Phenomena*, Springer Lecture Notes in
  Physics **186** (1983). Upper critical dimension and dangerous
  irrelevant variables; breakdown of hyperscaling at ``d > 4``.
- F. Kos, D. Poland, D. Simmons-Duffin, *Precision islands in the
  Ising and O(N) models*, Phys. Rev. D **93**, 036404 (2016).
  Precise 3D Ising exponents from the conformal bootstrap; used
  for the 3D check in Step 8.
- J. Cardy, *Scaling and Renormalization in Statistical Physics*,
  Cambridge University Press (1996), Ch. 3. Textbook treatment of
  scaling relations from the RG point of view.

## Used by

- [Ising universality class](../universalities/ising.md) — the
  exponent table is stored in `src/universalities/Ising2D.jl` using
  `Rational{Int}` values that satisfy all four scaling relations
  exactly; similarly for the 3D and mean-field entries.
- [`ising-cft-primary-operators`](ising-cft-primary-operators.md) —
  provides the two CFT scaling dimensions $(\Delta_{\sigma},
  \Delta_{\varepsilon}) = (1/8, 1)$ from which all six
  thermodynamic exponents are derived in Step 6 above.
- [`yang-magnetization-toeplitz`](yang-magnetization-toeplitz.md) —
  microscopic derivation of ``\beta = 1/8`` from the Toeplitz
  determinant, cross-checked against ``\beta = \nu\,\Delta_{\sigma}``
  here.
- [Percolation universality](../universalities/percolation.md),
  [Potts universality](../universalities/potts.md),
  [KPZ universality](../universalities/kpz.md), [O(N)
  universality](../universalities/on-models.md),
  [Mean-field universality](../universalities/mean-field.md) —
  same four scaling relations, different exponent tuples; QAtlas
  verifies each with `Rational{Int}` arithmetic.
