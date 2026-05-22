# Lieb Lattice Flat Band at ``E = 0``

## Main result

For the nearest-neighbour tight-binding model on the Lieb lattice
(three sublattices ``A, B, C`` per unit cell, ``A`` at square corners,
``B`` and ``C`` at the two edge midpoints) with hopping amplitude ``t``,
the ``3\times 3`` Bloch Hamiltonian has one **perfectly flat**
eigenvalue across the entire Brillouin zone,

```math
\boxed{\;
E_{\rm flat}(\mathbf{k}) \;=\; 0
\qquad\text{for every }\mathbf{k} \in \mathrm{BZ},
\;}
```

and two dispersive bands symmetric about ``E = 0``,

```math
\boxed{\;
E_{\pm}(\mathbf{k}) \;=\; \pm\,2\,t\,
      \sqrt{\,\cos^{2}\!\tfrac{\theta_{1}}{2}
             + \cos^{2}\!\tfrac{\theta_{2}}{2}\,},
\qquad \theta_{j} \equiv \mathbf{k}\cdot\mathbf{a}_{j}.
\;}
```

The dispersive bands touch the flat band at a single point — the
**M point** ``\mathbf{k} = (\pi/a, \pi/a)`` at the Brillouin-zone
corner — where ``\cos(\theta_{1}/2) = \cos(\theta_{2}/2) = 0`` and
``E_{\pm}(M) = 0``. Expansion near ``M`` yields a **linear (Dirac-
like)** touching of the three bands, all meeting at ``E = 0``.

The flat band's origin is the **bipartite sublattice imbalance** of
the Lieb lattice: the ``A`` sublattice has ``1`` site per unit cell
while ``B\cup C`` has ``2``. Lieb's theorem (1989) guarantees at least
``|A| - |B\cup C|`` zero-energy eigenvectors in magnitude, which for
``N_{c}`` unit cells is exactly ``N_{c}`` — one per ``\mathbf{k}``, i.e.
a full flat band.

---

## Setup

### Lattice geometry

The Lieb lattice is a decoration of the square lattice: starting
from a plain square of primitive vectors ``\mathbf{a}_{1} = a\,(1, 0)``,
``\mathbf{a}_{2} = a\,(0, 1)``, add one site at the midpoint of each
horizontal edge and one at the midpoint of each vertical edge.
Denote these

- ``A`` (corner / "copper"): position ``\mathbf{r}_{A} = (0, 0)``,
- ``B`` (horizontal-edge midpoint / "oxygen ``x``"): $\mathbf{r}_{B}
  = \mathbf{a}_{1}/2 = (a/2, 0)$,
- ``C`` (vertical-edge midpoint / "oxygen ``y``"): $\mathbf{r}_{C}
  = \mathbf{a}_{2}/2 = (0, a/2)$.

(The naming anticipates the CuO``_{2}`` plane of cuprates, where ``A``
carries the Cu orbital and ``B, C`` carry the O``_{x}, \text{O}_{y}``
orbitals.) Set ``a = 2`` so that the NN bond length $|\mathbf{r}_{B}
- \mathbf{r}_{A}| = 1``; then ``\mathbf{a}_{1} = (2, 0)$,
``\mathbf{a}_{2} = (0, 2)``.

Each ``A``-site at ``\mathbf{R}`` has four NN: two ``B``-sites along
``\pm \hat{x}`` and two ``C``-sites along ``\pm \hat{y}``. Each ``B``-site
has only two NN: the two ``A``-sites flanking it along ``\hat{x}``.
Similarly each ``C``-site has two NN (the two ``A``-sites along
``\hat{y}``). There is **no** direct ``B``–``C`` hopping — ``B`` and ``C``
sites are **not** nearest neighbours (the ``B``–``C`` distance is
``\sqrt{(a/2)^{2} + (a/2)^{2}} = a/\sqrt{2}``, which is a second-
neighbour distance).

Explicitly, from the ``A``-site at ``\mathbf{R}``:

```math
A_{\mathbf{R}} \leftrightarrow
\begin{cases}
B_{\mathbf{R}} & \text{(same cell)}\\
B_{\mathbf{R} - \mathbf{a}_{1}} & \text{(cell } -\mathbf{a}_{1}\text{)} \\
C_{\mathbf{R}} & \text{(same cell)}\\
C_{\mathbf{R} - \mathbf{a}_{2}} & \text{(cell } -\mathbf{a}_{2}\text{)}
\end{cases}
\tag{1}
```

The lattice is **bipartite** with partition ``(A) \sqcup (B\cup C)``:
every NN bond connects an ``A``-site to a ``B``- or ``C``-site.

### Bipartite structure and sublattice imbalance

On an ``L_{1}\times L_{2}`` unit-cell torus (``N_{c} = L_{1}L_{2}``),

```math
|A| = N_{c},\qquad |B| = N_{c},\qquad |C| = N_{c},
\qquad |B\cup C| = 2 N_{c}.
```

The imbalance ``|B\cup C| - |A| = N_{c}`` is the key input for Lieb's
theorem (Step 5), which predicts at least ``N_{c}`` zero-energy
eigenvectors — the entire flat band.

### Hamiltonian

```math
H \;=\; -\,t\sum_{\langle i, j\rangle}
   \bigl(c^{\dagger}_{i}c_{j} + c^{\dagger}_{j}c_{i}\bigr),
```

with the sum over NN ``A``–``B`` and ``A``–``C`` bonds only (no ``B``–``C``).

### Goal

Derive the ``3\times 3`` Bloch Hamiltonian, find the exact flat-band
eigenvector, and present the bipartite-imbalance argument that
explains its origin.

---

## Calculation

### Step 1 — Fourier transform to the Bloch Hamiltonian

Use the site-centered Fourier convention (as in the
[kagome note](bloch-kagome-flat-band.md), Step 2):

```math
c_{\alpha, \mathbf{k}}
 \;=\; \frac{1}{\sqrt{N_{c}}}\sum_{\mathbf{R}}
         e^{-i\mathbf{k}\cdot(\mathbf{R} + \mathbf{r}_{\alpha})}
         c_{\alpha, \mathbf{R}}.
```

Each NN bond between ``\alpha`` at ``\mathbf{r}_{\alpha} + \mathbf{R}``
and ``\beta`` at ``\mathbf{r}_{\beta} + \mathbf{R}'`` contributes
``-t\,e^{-i\mathbf{k}\cdot\mathbf{d}}`` in the Bloch Hamiltonian
``(\alpha, \beta)`` block, with $\mathbf{d} \equiv (\mathbf{r}_{\beta}
+ \mathbf{R}') - (\mathbf{r}_{\alpha} + \mathbf{R})$ the bond
vector. A Hermitian-conjugate pair sums to $-2t\cos(\mathbf{k}\cdot
\mathbf{d})$.

For the two ``A``–``B`` bonds of (1): same-cell $\mathbf{d} =
+\mathbf{a}_{1}/2``, cross-cell ``\mathbf{d} = -\mathbf{a}_{1}/2$.
Summing,

```math
\widetilde{H}_{AB}(\mathbf{k})
 \;=\; -\,t\,\bigl(e^{-i\mathbf{k}\cdot\mathbf{a}_{1}/2}
           + e^{+i\mathbf{k}\cdot\mathbf{a}_{1}/2}\bigr)
 \;=\; -\,2 t\cos(\theta_{1}/2).
```

Similarly ``\widetilde{H}_{AC}(\mathbf{k}) = -2t\cos(\theta_{2}/2)``.
The ``(B, C)`` entry is zero (no direct ``B``–``C`` hopping). Assembling,

```math
\boxed{\;
\widetilde{H}(\mathbf{k})
 \;=\; -\,2 t\,\begin{pmatrix}
 0 & \cos(\theta_{1}/2) & \cos(\theta_{2}/2) \\
 \cos(\theta_{1}/2) & 0 & 0 \\
 \cos(\theta_{2}/2) & 0 & 0
\end{pmatrix}.
\;}
\tag{2}
```

The matrix is real symmetric and has the generic **bipartite block
structure**

```math
\widetilde{H}(\mathbf{k})
 \;=\; \begin{pmatrix} 0 & \mathbf{q}^{T}(\mathbf{k}) \\
   \mathbf{q}(\mathbf{k}) & \mathbf{0}_{2\times 2}\end{pmatrix},
\qquad
\mathbf{q}(\mathbf{k}) \;=\; -2t\,\begin{pmatrix}
   \cos(\theta_{1}/2) \\ \cos(\theta_{2}/2)\end{pmatrix}
\in \mathbb{R}^{2}.
\tag{3}
```

The ``2\times 2`` lower-right zero block reflects the absence of
direct ``B``–``C`` hopping.

### Step 2 — Secular equation

Set ``\lambda = -2t\,\mu`` so that ``\mu`` is an eigenvalue of
``M\equiv\widetilde{H}/(-2t)``. By cofactor expansion along the first
row of ``M - \mu\mathbb{I}`` with ``s_{1} = \cos(\theta_{1}/2)``,
``s_{2} = \cos(\theta_{2}/2)``,

```math
\det(M - \mu\mathbb{I})
 = -\mu\,\det\!\begin{pmatrix}-\mu & 0 \\ 0 & -\mu\end{pmatrix}
   - s_{1}\,\det\!\begin{pmatrix}s_{1} & 0 \\ s_{2} & -\mu\end{pmatrix}
   + s_{2}\,\det\!\begin{pmatrix}s_{1} & -\mu \\ s_{2} & 0\end{pmatrix}
```
```math
= -\mu^{3} + \mu\,s_{1}^{2} + \mu\,s_{2}^{2}
 = -\mu\,\bigl[\,\mu^{2} - (s_{1}^{2} + s_{2}^{2})\,\bigr].
```

Setting the determinant to zero,

```math
\boxed{\;
\mu\,\bigl[\,\mu^{2} - (s_{1}^{2} + s_{2}^{2})\,\bigr] \;=\; 0.
\;}
\tag{4}
```

The three roots are

```math
\mu_{0} = 0,
\qquad
\mu_{\pm} = \pm\,\sqrt{s_{1}^{2} + s_{2}^{2}}.
```

Translating back to ``\widetilde{H}`` eigenvalues ``\lambda = -2t\mu``,

```math
E_{\rm flat} \;=\; 0
\qquad\text{(all }\mathbf{k}\text{)},
\qquad
E_{\pm}(\mathbf{k}) \;=\; \mp\,2 t\,\sqrt{s_{1}^{2} + s_{2}^{2}}
 \;=\; \pm\,2 t\,\sqrt{\cos^{2}\!\tfrac{\theta_{1}}{2}
               + \cos^{2}\!\tfrac{\theta_{2}}{2}}.
\tag{5}
```

The flat band sits at ``E = 0`` throughout the BZ, independent of
``\mathbf{k}``.

### Step 3 — Explicit flat-band eigenvector

Solve ``\widetilde{H}(\mathbf{k})\,\mathbf{v} = 0``. In components
``\mathbf{v} = (v_{A}, v_{B}, v_{C})^{T}``, (2) gives

```math
\text{row }A: \quad -2t\bigl[s_{1} v_{B} + s_{2} v_{C}\bigr] = 0,
```
```math
\text{row }B: \quad -2t\,s_{1} v_{A} = 0,
```
```math
\text{row }C: \quad -2t\,s_{2} v_{A} = 0.
```

For generic ``\mathbf{k}`` with ``s_{1}, s_{2} \ne 0``, rows ``B`` and ``C``
force ``v_{A} = 0``. Row ``A`` then requires $s_{1} v_{B} + s_{2} v_{C}
= 0``, a single linear constraint on ``(v_{B}, v_{C})$. A convenient
normalisation picks

```math
\boxed{\;
|\psi_{\rm flat}(\mathbf{k})\rangle \;=\;
 \frac{1}{\sqrt{s_{1}^{2} + s_{2}^{2}}}\,
 \bigl(\,0,\;\; -s_{2},\;\; +s_{1}\,\bigr)^{T},
\;}
\tag{6}
```

so that ``\langle\psi_{\rm flat}|\psi_{\rm flat}\rangle = 1`` for all
``\mathbf{k}`` except the single ``M``-point $\theta_{1} = \theta_{2} =
\pi`` where both ``s_{j} = 0$ (see Step 6). The flat-band
wavefunction **has zero weight on the ``A`` sublattice for every
``\mathbf{k}``** — this is the algebraic expression of the bipartite
imbalance argument (Step 5).

### Step 4 — Dispersive bands and the ``E_{\pm}`` eigenvectors

The ``E_{\pm}`` eigenvectors satisfy $\widetilde{H}\,\mathbf{v}_{\pm}
= \pm E \mathbf{v}_{\pm}`` with ``E = 2t\sqrt{s_{1}^{2} + s_{2}^{2}}$.
From rows ``B`` and ``C``,

```math
-2t\,s_{1}\,v_{A} = \pm E\,v_{B},
\qquad
-2t\,s_{2}\,v_{A} = \pm E\,v_{C},
```

giving ``v_{B} = -2t s_{1} v_{A}/(\pm E)``, $v_{C} = -2t s_{2}
v_{A}/(\pm E)$. Normalising,

```math
|\psi_{\pm}(\mathbf{k})\rangle \;=\;
 \frac{1}{\sqrt{2(s_{1}^{2} + s_{2}^{2})}}\,
 \begin{pmatrix}\pm\sqrt{s_{1}^{2} + s_{2}^{2}} \\
                -s_{1} \\ -s_{2}\end{pmatrix}.
\tag{7}
```

Both ``\psi_{\pm}`` have non-zero weight on the ``A`` sublattice; only
the flat-band state (6) is ``A``-null. The three bands are orthogonal
for every ``\mathbf{k}`` where ``E(\mathbf{k}) > 0``.

### Step 5 — Bipartite-imbalance argument (Lieb 1989)

The flat band of (5) is no accident. It follows from a general
theorem on bipartite lattices: if a tight-binding Hamiltonian has
the block structure ``M = \begin{pmatrix}0 & Q\\ Q^{T} & 0\end{pmatrix}``
with ``Q \in \mathbb{R}^{N_{1}\times N_{2}}`` (``N_{1}`` sites in one
sublattice, ``N_{2}`` in the other), then the **at-least-$|N_{1} -
N_{2}|``-degeneracy** at ``E = 0$ follows from

```math
M\,\mathbf{v} = 0
\quad\Longleftrightarrow\quad
Q\,\mathbf{v}_{2} = 0\quad\text{and}\quad Q^{T}\,\mathbf{v}_{1} = 0.
```

By rank-nullity, $\dim\ker Q = N_{2} - \mathrm{rank}(Q) \ge
N_{2} - \min(N_{1}, N_{2}) = \max(0, N_{2} - N_{1})$, so there are
at least ``\max(0, N_{2} - N_{1})`` zero-energy eigenvectors
supported on sublattice 2. Similarly at least
``\max(0, N_{1} - N_{2})`` such vectors on sublattice 1. In total, at
least ``|N_{1} - N_{2}|`` zero-energy eigenvectors regardless of the
detailed structure of ``Q``.

For the Lieb lattice, sublattice 1 is ``A`` (``N_{1} = N_{c}``) and
sublattice 2 is ``B\cup C`` (``N_{2} = 2 N_{c}``). The imbalance is
``|N_{2} - N_{1}| = N_{c}``, giving ``\ge N_{c}`` zero-energy states on
``B\cup C`` (consistent with our real-space expectation
"``v_{A} = 0``"). Distributing these over the ``N_{c}`` values of
``\mathbf{k}`` gives one zero mode per momentum, i.e. **a complete
flat band**. This is Lieb's theorem (Phys. Rev. Lett. **62**, 1201
(1989)) specialised to the Lieb-lattice adjacency matrix.

The argument generalises: any bipartite lattice with unequal
sublattice sizes has at least ``|N_{1} - N_{2}|`` zero-energy modes,
and if the lattice is connected (every site reachable via NN
bonds) the degeneracy saturates the bound. For the Lieb lattice at
half-filling this macroscopic degeneracy produces **flat-band
ferromagnetism** in the Hubbard model (Lieb 1989, Mielke 1991,
Tasaki 1998).

### Step 6 — Band-touching at the ``M``-point

The flat band and the dispersive bands (5) coincide at ``E = 0``
wherever ``s_{1}^{2} + s_{2}^{2} = 0``, i.e. $\cos(\theta_{1}/2) =
\cos(\theta_{2}/2) = 0``, i.e. ``\theta_{1} = \theta_{2} = \pi$. This
is the BZ corner (the M-point) at ``\mathbf{k} = (\pi/a, \pi/a)``.

Expand (5) near ``M``. With our normalisation $\mathbf{a}_{j} = 2
\hat{e}_{j}``, we have ``\theta_{j} = \mathbf{k}\cdot\mathbf{a}_{j} =
2 k_{j}``; the ``M``-point in the ``\mathbf{k}$-variable is
``\mathbf{M} = (\pi/2, \pi/2)``. Set $\mathbf{q} = \mathbf{k} -
\mathbf{M}``, so ``\theta_{j} = \pi + 2 q_{j}$ and

```math
\cos(\theta_{j}/2) \;=\; \cos(\pi/2 + q_{j})
 \;=\; -\sin q_{j} \;\approx\; -q_{j} + O(q_{j}^{3}).
```

Therefore

```math
s_{1}^{2} + s_{2}^{2} \;\approx\; q_{1}^{2} + q_{2}^{2} \;=\; |\mathbf{q}|^{2},
```

and

```math
E_{\pm}(\mathbf{M} + \mathbf{q})
 \;=\; \pm\,2t\,\sqrt{q_{1}^{2} + q_{2}^{2}}
 \;=\; \pm\,2t\,|\mathbf{q}| \;+\; O(|\mathbf{q}|^{3}),
\tag{8}
```

so the dispersive bands touch the flat band **linearly** — a
three-way Dirac-like touching at ``M``. The two dispersive bands plus
the flat band all meet at ``E = 0``, with the dispersive cones
meeting the plane tangentially along a cone and a flat pancake
respectively. The Fermi velocity is ``v_{F} = 2t``.

### Step 7 — Limiting-case and finite-size checks

**(i) Half-filling.** With ``3 N_{c}`` sites, half-filling is
``3 N_{c}/2`` electrons. If ``N_{c}`` is odd the count is half-integer
so one assumes a grand-canonical ensemble. At ``T = 0`` the ``N_{c}``
states of the lower dispersive band (``E_{-} < 0``) are filled; the
flat band at ``E = 0`` is half-filled (in the spinless case) or
fully-filled per spin (in the spinful case, with additional flat
spin degeneracy). This macroscopic degeneracy is the input for
flat-band ferromagnetism (Lieb 1989).

**(ii) ``\Gamma``-point.** At ``\mathbf{k} = (0, 0)``: $s_{1} = s_{2}
= 1``, ``s_{1}^{2} + s_{2}^{2} = 2``, so ``E_{\pm}(\Gamma) = \pm 2t
\sqrt{2}``. Three-band spectrum at ``\Gamma``: ``\{-2t\sqrt{2},\, 0,\,
+2t\sqrt{2}\}``, all distinct (bandwidth ``= 4t\sqrt{2}$).

**(iii) Bipartite chiral symmetry.** The block form (3) makes the
spectrum symmetric about ``E = 0``: for every dispersive eigenvalue
``+E(\mathbf{k})`` there is ``-E(\mathbf{k})``. The sublattice operator
``\Sigma = \mathrm{diag}(+1, -1, -1)`` (positive on ``A``,
negative on ``B\cup C``) satisfies $\Sigma \widetilde{H}(\mathbf{k})
\Sigma = -\widetilde{H}(\mathbf{k})$, the algebraic expression of
bipartiteness.

**(iv) Honeycomb vs Lieb comparison.** The honeycomb lattice (see
[bloch-honeycomb-dispersion](bloch-honeycomb-dispersion.md)) is also
bipartite but with equal sublattice sizes ``|A| = |B|``, so Lieb's
theorem gives zero zero-mode lower bound — and indeed the honeycomb
spectrum has Dirac touching but no flat band. The Lieb lattice
differs precisely by its ``1 : 2`` sublattice imbalance.

**(v) Kagome vs Lieb comparison.** The [kagome flat band](bloch-kagome-flat-band.md)
sits at ``+2t``, not at ``E = 0``, and arises from destructive
interference on hexagonal plaquettes — not from a bipartite
sublattice imbalance (the kagome lattice is **not** bipartite).
The two flat-band mechanisms are genuinely distinct paradigms.

---

## References

- E. H. Lieb, *Two theorems on the Hubbard model*, Phys. Rev.
  Lett. **62**, 1201 (1989). Bipartite sublattice imbalance →
  ground-state spin ``S = |N_{A} - N_{B}|/2``; flat-band
  ferromagnetism in the half-filled Hubbard model.
- A. Mielke, *Ferromagnetism in the Hubbard model on line graphs
  and further considerations*, J. Phys. A **24**, 3311 (1991).
  Flat bands on line graphs; Lieb-lattice case as a decorated
  bipartite lattice.
- H. Tasaki, *From Nagaoka's ferromagnetism to flat-band
  ferromagnetism and beyond*, Prog. Theor. Phys. **99**, 489 (1998).
  Rigorous review of flat-band magnetism mechanisms including the
  Lieb / Mielke constructions.
- Y.-F. Wang, C.-D. Gong, *Frustration induced non-trivial flat
  bands in Lieb lattice*, Phys. Rev. B **85**, 214428 (2012).
  Modern treatment including perturbations.

## Used by

- [Lieb tight-binding model](../models/quantum/tightbinding/lieb.md) —
  `fetch(QAtlas.Lieb(; t, Lx, Ly), TightBindingSpectrum(), …)`
  returns the finite-size three-band spectrum, with exact zeros at
  ``N_{c}`` momenta (the full flat band) and the dispersive ``\pm E(\mathbf{k})``.
- [Honeycomb dispersion note](bloch-honeycomb-dispersion.md) —
  bipartite lattice with equal sublattice sizes gives a linear
  Dirac cone but no flat band; the Lieb ``1:2`` imbalance produces
  the flat band.
- [Kagome flat-band note](bloch-kagome-flat-band.md) — a second,
  fundamentally different flat-band paradigm: destructive
  interference on a non-bipartite three-sublattice lattice.
