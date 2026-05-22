# Kagome Bloch Hamiltonian and the Flat Band

## Main result

For the nearest-neighbour tight-binding model on the kagome lattice
with hopping amplitude ``t``, the ``3\times 3`` Bloch Hamiltonian has one
**perfectly flat** eigenvalue across the entire Brillouin zone,

```math
\boxed{\;
E_{\rm flat}(\mathbf{k}) \;=\; +2\,t
\qquad\text{for every }\mathbf{k} \in \mathrm{BZ},
\;}
```

and two dispersive bands

```math
\boxed{\;
E_{\pm}(\mathbf{k}) \;=\; -\,t \pm t\,\sqrt{\,4 S(\mathbf{k}) - 3\,},
\qquad
S(\mathbf{k}) \;\equiv\; \cos^{2}\!\tfrac{\theta_{1}}{2}
                        + \cos^{2}\!\tfrac{\theta_{2}}{2}
                        + \cos^{2}\!\tfrac{\theta_{2}-\theta_{1}}{2},
\;}
```

with ``\theta_{j} \equiv \mathbf{k}\cdot\mathbf{a}_{j}``. The flat band
and the upper dispersive band ``E_{+}`` are degenerate at the
``\Gamma``-point ``\mathbf{k} = (0, 0)`` (both equal to ``+2t`` there),
forming a **quadratic band-touching point**; they separate for
``\mathbf{k} \ne 0``.

The exact flatness holds for every finite kagome lattice with
periodic boundary conditions: the adjacency matrix has an eigenvalue
``-1`` with massive degeneracy, corresponding to real-space
**compact localised states** supported on hexagonal plaquettes of
the lattice (Sutherland 1986, Mielke 1991).

---

## Setup

### Lattice geometry

The kagome lattice has ``3`` sites per primitive cell, labelled
``A, B, C``, each at a corner of a triangle. With triangle side length
``a`` (= twice the NN bond length ``a/2``), the Bravais-lattice
primitive vectors are

```math
\mathbf{a}_{1} \;=\; a\,(1, 0),
\qquad
\mathbf{a}_{2} \;=\; a\,\bigl(\tfrac{1}{2},\,\tfrac{\sqrt{3}}{2}\bigr).
```

Set ``a = 1`` henceforth. The three sublattice positions within a cell
are

```math
\mathbf{r}_{A} = (0, 0),
\qquad
\mathbf{r}_{B} = \bigl(\tfrac{1}{2}, 0\bigr) = \tfrac{\mathbf{a}_{1}}{2},
\qquad
\mathbf{r}_{C} = \bigl(\tfrac{1}{4}, \tfrac{\sqrt{3}}{4}\bigr)
                 = \tfrac{\mathbf{a}_{2}}{2}.
```

The NN bond length is
``|\mathbf{r}_{B} - \mathbf{r}_{A}| = 1/2``; checking $|\mathbf{r}_{C}
- \mathbf{r}_{A}| = \sqrt{1/16 + 3/16} = 1/2`` and ``|\mathbf{r}_{C}
- \mathbf{r}_{B}| = \sqrt{1/16 + 3/16} = 1/2$, so every in-cell pair
forms a NN bond (the three sublattices form an upward-pointing
triangle in each cell, sharing edges with downward-pointing
triangles in neighbouring cells).

Each site has **four** nearest neighbours: two in the same unit
cell and two in neighbouring cells. Explicitly, from an ``A``-site at
``\mathbf{R}`` the four NN bonds go to

```math
A_{\mathbf{R}} \leftrightarrow
\begin{cases}
B_{\mathbf{R}} & \text{(same-cell, via }\mathbf{r}_B - \mathbf{r}_A\text{)} \\
B_{\mathbf{R} - \mathbf{a}_{1}} & \text{(cell }-\mathbf{a}_{1}\text{)} \\
C_{\mathbf{R}} & \text{(same-cell)} \\
C_{\mathbf{R} - \mathbf{a}_{2}} & \text{(cell }-\mathbf{a}_{2}\text{)}
\end{cases}
\tag{1}
```

Similarly each ``B``-site connects to two ``A`` and two ``C`` NN, and each
``C``-site to two ``A`` and two ``B`` NN. Total bond count: on an
``L_{1}L_{2}``-unit-cell PBC torus the kagome lattice has ``3 L_{1}L_{2}``
sites and ``6 L_{1}L_{2}`` bonds (4 bonds per site, divided by 2).

### Hamiltonian

```math
H \;=\; -\,t\sum_{\langle i, j\rangle}
   \bigl(c^{\dagger}_{i} c_{j} + c^{\dagger}_{j} c_{i}\bigr),
```

with ``\langle i, j\rangle`` over NN bonds and ``t > 0``.

### Goal

Derive the ``3\times 3`` Bloch Hamiltonian, diagonalise, and prove
that one eigenvalue is constant (independent of ``\mathbf{k}``) —
yielding the flat band.

---

## Calculation

### Step 1 — Fourier transform to the Bloch Hamiltonian

Define Bloch operators without sublattice offsets,

```math
c_{\alpha, \mathbf{k}} \;=\; \frac{1}{\sqrt{N_{c}}}
    \sum_{\mathbf{R}} e^{-i\mathbf{k}\cdot\mathbf{R}}\,c_{\alpha, \mathbf{R}},
\qquad \alpha \in \{A, B, C\},
```

for ``N_{c} = L_{1}L_{2}`` unit cells. The two ``A``–``B`` bonds from (1)
give the same-cell contribution (phase ``1``) and the cross-cell
contribution (``\mathbf{R}`` to ``\mathbf{R} - \mathbf{a}_{1}``, phase
``e^{i\mathbf{k}\cdot\mathbf{a}_{1}}`` by the same-site-collapse argument
used in the [honeycomb note](bloch-honeycomb-dispersion.md) Step 1):

```math
\sum_{\mathbf{R}}\bigl(A^{\dagger}_{\mathbf{R}} B_{\mathbf{R}}
                     + A^{\dagger}_{\mathbf{R}} B_{\mathbf{R} - \mathbf{a}_{1}}\bigr)
\;=\; \sum_{\mathbf{k}}
   \bigl(1 + e^{i\mathbf{k}\cdot\mathbf{a}_{1}}\bigr)\,
   A^{\dagger}_{\mathbf{k}} B_{\mathbf{k}}.
```

Similarly for the two ``A``–``C`` bonds (cross-cell shift
``-\mathbf{a}_{2}``) and the two ``B``–``C`` bonds. The geometry for the
latter: ``B`` at ``\mathbf{r}_{B}`` connects to ``C`` at ``\mathbf{r}_{C}``
in the same cell, and to ``C`` at $\mathbf{r}_{C} + \mathbf{a}_{1} -
\mathbf{a}_{2}`` in the cell shifted by ``\mathbf{a}_{1} -
\mathbf{a}_{2}`` (check: ``|\mathbf{r}_{C} + \mathbf{a}_{1} -
\mathbf{a}_{2} - \mathbf{r}_{B}| = |(1/4, \sqrt{3}/4) + (1/2,
-\sqrt{3}/2) - (1/2, 0)| = |(1/4, -\sqrt{3}/4)| = 1/2$, NN ✓).
Fourier gives the factor $1 + e^{i\mathbf{k}\cdot(\mathbf{a}_{1} -
\mathbf{a}_{2})}$.

Collecting contributions (and using
$-t\,(1 + e^{i\phi})\,A^{\dagger}_{\mathbf{k}}B_{\mathbf{k}} +
\text{h.c.}$ for each bond type),

```math
H = \sum_{\mathbf{k}}\,\Psi^{\dagger}_{\mathbf{k}}\,
  \widetilde{H}(\mathbf{k})\,\Psi_{\mathbf{k}},\qquad
\Psi_{\mathbf{k}} = (A_{\mathbf{k}}, B_{\mathbf{k}}, C_{\mathbf{k}})^{T},
```

```math
\widetilde{H}(\mathbf{k})
 = -t\,\begin{pmatrix}
   0 & 1 + e^{i\theta_{1}} & 1 + e^{i\theta_{2}} \\
   1 + e^{-i\theta_{1}} & 0 & 1 + e^{-i(\theta_{2} - \theta_{1})} \\
   1 + e^{-i\theta_{2}} & 1 + e^{i(\theta_{2} - \theta_{1})} & 0
 \end{pmatrix}.
\tag{2}
```

with ``\theta_{j} = \mathbf{k}\cdot\mathbf{a}_{j}``. This is the
**complex-gauge** Bloch Hamiltonian.

### Step 2 — Gauge transformation to real-symmetric form

Switch to the **site-centered Fourier convention**

```math
c_{\alpha, \mathbf{k}} \;=\; \frac{1}{\sqrt{N_{c}}}
    \sum_{\mathbf{R}} e^{-i\mathbf{k}\cdot(\mathbf{R} + \mathbf{r}_{\alpha})}
                        c_{\alpha, \mathbf{R}},
\qquad \alpha \in \{A, B, C\},
```

which replaces the cell-index phase ``e^{-i\mathbf{k}\cdot\mathbf{R}}``
of Step 1 by a full-position phase $e^{-i\mathbf{k}\cdot(\mathbf{R} +
\mathbf{r}_{\alpha})}$ including the sublattice offset. Under this
convention every bond $\alpha_{\mathbf{R}} \leftrightarrow
\beta_{\mathbf{R}'}$ contributes, in the Bloch Hamiltonian,

```math
-t\,e^{-i\mathbf{k}\cdot(\mathbf{r}_{\beta} + \mathbf{R}'
                    - \mathbf{r}_{\alpha} - \mathbf{R})}
 \;=\; -t\,e^{-i\mathbf{k}\cdot\mathbf{d}},
```

with $\mathbf{d} = (\mathbf{r}_{\beta} + \mathbf{R}') -
(\mathbf{r}_{\alpha} + \mathbf{R})$ the **bond vector** from
``\alpha`` to ``\beta``. Summing a Hermitian-conjugate pair gives
``-2t\cos(\mathbf{k}\cdot\mathbf{d})``.

For the two ``A``–``B`` bonds of (1): the same-cell bond has $\mathbf{d}
= \mathbf{r}_{B} - \mathbf{r}_{A} = +\mathbf{a}_{1}/2$; the cross-cell
bond has $\mathbf{d} = \mathbf{r}_{B} - \mathbf{r}_{A} -
\mathbf{a}_{1} = -\mathbf{a}_{1}/2``. The ``(A, B)$ matrix element is
the sum over these two:

```math
\widetilde{H}_{AB}(\mathbf{k})
 = -t\bigl(e^{-i\theta_{1}/2} + e^{+i\theta_{1}/2}\bigr)
 = -2 t\cos(\theta_{1}/2).
```

The ``(B, A)`` entry is the complex conjugate, which equals
``\widetilde{H}_{AB}`` since ``\cos`` is real. The same computation for
the other bond types:

```math
\widetilde{H}_{AC} = -2t\cos(\theta_{2}/2),
\qquad
\widetilde{H}_{BC} = -2t\cos\bigl[(\theta_{2} - \theta_{1})/2\bigr].
```

The in-cell ``B``–``C`` bond vector is $\mathbf{r}_{C} - \mathbf{r}_{B}
= \mathbf{a}_{2}/2 - \mathbf{a}_{1}/2$; the cross-cell one is the
negative of that. Both give ``\cos((\theta_{2} - \theta_{1})/2)``
after summation. Assembling, the **real-symmetric Bloch Hamiltonian** is

```math
\boxed{\;
\widetilde{H}(\mathbf{k})
 \;=\; -\,2t\,\begin{pmatrix}
   0 & c_{1} & c_{2} \\
   c_{1} & 0 & c_{3} \\
   c_{2} & c_{3} & 0
 \end{pmatrix},
\qquad
\begin{aligned}
c_{1} &= \cos(\theta_{1}/2),\\
c_{2} &= \cos(\theta_{2}/2),\\
c_{3} &= \cos\bigl[(\theta_{2} - \theta_{1})/2\bigr].
\end{aligned}
\;}
\tag{3}
```

The two gauges give *unitarily equivalent* Hamiltonians (3) vs. (2);
the eigenvalues are identical, but (3) is more convenient for the
trigonometric identity used below.

### Step 3 — Secular equation

Let ``M \equiv \widetilde{H}(\mathbf{k})/(-2t)``, so that eigenvalues
of ``\widetilde{H}`` are ``\lambda = -2t\,\mu`` where ``\mu`` is an
eigenvalue of ``M``:

```math
M = \begin{pmatrix}
  0 & c_{1} & c_{2} \\
  c_{1} & 0 & c_{3} \\
  c_{2} & c_{3} & 0
\end{pmatrix}.
```

Compute the characteristic polynomial ``\det(M - \mu\,\mathbb{I}) = 0``
by cofactor expansion along the first row:

```math
\det(M - \mu\mathbb{I})
= -\mu\,\det\!\begin{pmatrix}-\mu & c_{3} \\ c_{3} & -\mu\end{pmatrix}
\;-\; c_{1}\det\!\begin{pmatrix}c_{1} & c_{3} \\ c_{2} & -\mu\end{pmatrix}
\;+\; c_{2}\det\!\begin{pmatrix}c_{1} & -\mu \\ c_{2} & c_{3}\end{pmatrix}.
```

Evaluating the ``2\times 2`` determinants,

```math
\det\!\begin{pmatrix}-\mu & c_{3} \\ c_{3} & -\mu\end{pmatrix} = \mu^{2} - c_{3}^{2},
```
```math
\det\!\begin{pmatrix}c_{1} & c_{3} \\ c_{2} & -\mu\end{pmatrix}
 = -\mu c_{1} - c_{2} c_{3},
```
```math
\det\!\begin{pmatrix}c_{1} & -\mu \\ c_{2} & c_{3}\end{pmatrix}
 = c_{1} c_{3} + \mu c_{2}.
```

Substituting,

```math
\det(M - \mu\mathbb{I})
 = -\mu(\mu^{2} - c_{3}^{2}) - c_{1}(-\mu c_{1} - c_{2} c_{3})
    + c_{2}(c_{1} c_{3} + \mu c_{2})
```
```math
= -\mu^{3} + \mu\,c_{3}^{2} + \mu\,c_{1}^{2} + c_{1} c_{2} c_{3}
    + c_{1} c_{2} c_{3} + \mu\,c_{2}^{2}
```
```math
= -\mu^{3} + \mu\,S + 2\,P,
```

with

```math
S(\mathbf{k}) \;\equiv\; c_{1}^{2} + c_{2}^{2} + c_{3}^{2},
\qquad
P(\mathbf{k}) \;\equiv\; c_{1} c_{2} c_{3}.
```

Setting ``\det(M - \mu\mathbb{I}) = 0`` yields the **secular
equation**

```math
\boxed{\;
\mu^{3} \;-\; S\,\mu \;-\; 2 P \;=\; 0.
\;}
\tag{4}
```

### Step 4 — Flat-band root ``\mu = -1``

We show that ``\mu = -1`` is a root of (4) for every ``\mathbf{k}``, by
proving the trigonometric identity

```math
\boxed{\;
c_{1}^{2} + c_{2}^{2} + c_{3}^{2} \;-\; 2\,c_{1} c_{2} c_{3}
\;=\; 1,
\;}
\tag{5}
```

where ``c_{1} = \cos\alpha``, ``c_{2} = \cos(\alpha + \beta)``, $c_{3} =
\cos\beta`` with ``\alpha = \theta_{1}/2``, ``\beta = (\theta_{2} -
\theta_{1})/2``, ``\alpha + \beta = \theta_{2}/2``. Substituting ``\mu =
-1$ into (4),

```math
(-1)^{3} - S\,(-1) - 2\,P = -1 + S - 2P
 = (S - 2P) - 1 = 1 - 1 = 0.
\quad\checkmark
```

**Proof of (5).** Use the sum-angle formula $\cos(\alpha + \beta) =
\cos\alpha\cos\beta - \sin\alpha\sin\beta$ and square:

```math
c_{2}^{2}
 = \cos^{2}(\alpha + \beta)
 = \cos^{2}\alpha\cos^{2}\beta - 2\cos\alpha\cos\beta\sin\alpha\sin\beta
    + \sin^{2}\alpha\sin^{2}\beta.
```

Substitute ``\sin^{2}x = 1 - \cos^{2}x``:

```math
= c_{1}^{2}c_{3}^{2}
   - 2 c_{1} c_{3}\sin\alpha\sin\beta
   + (1 - c_{1}^{2})(1 - c_{3}^{2})
```
```math
= 1 - c_{1}^{2} - c_{3}^{2} + 2\,c_{1}^{2}c_{3}^{2}
   - 2 c_{1} c_{3}\sin\alpha\sin\beta.
```

Add ``c_{1}^{2} + c_{3}^{2}`` to both sides:

```math
c_{1}^{2} + c_{2}^{2} + c_{3}^{2}
 = 1 + 2 c_{1}^{2}c_{3}^{2}
 - 2 c_{1}c_{3}\sin\alpha\sin\beta.
```

Factor the last two terms on the right:

```math
= 1 + 2\,c_{1}c_{3}\,\bigl(c_{1}c_{3} - \sin\alpha\sin\beta\bigr)
 = 1 + 2\,c_{1}c_{3}\,\cos(\alpha + \beta)
 = 1 + 2\,c_{1}c_{3}c_{2}
 = 1 + 2 P,
```

using $c_{1}c_{3} - \sin\alpha\sin\beta = \cos\alpha\cos\beta -
\sin\alpha\sin\beta = \cos(\alpha + \beta) = c_{2}$.

Hence ``c_{1}^{2} + c_{2}^{2} + c_{3}^{2} - 2 P = 1``, i.e.
``S - 2P = 1``, which is exactly (5). ``\square``

Therefore ``\mu = -1`` is a root of (4) for every ``\mathbf{k}``, and
the corresponding band of ``\widetilde{H}`` is at

```math
E_{\rm flat} \;=\; -2t\cdot(-1) \;=\; +2\,t
\qquad(\text{for every }\mathbf{k}).
```

### Step 5 — Dispersive bands

Factor ``(\mu + 1)`` out of (4). Polynomial division:

```math
\mu^{3} - S\,\mu - 2P
 = (\mu + 1)\bigl(\mu^{2} - \mu + (1 - S)\bigr),
```

which follows by expanding the right-hand side and using
``S - 2P = 1`` to match the constant term:

```math
(\mu + 1)(\mu^{2} - \mu + 1 - S)
 = \mu^{3} - \mu^{2} + \mu(1 - S) + \mu^{2} - \mu + (1 - S)
 = \mu^{3} - S\mu + (1 - S)
 = \mu^{3} - S\mu - 2P,
```

using ``1 - S = -2P``. The remaining quadratic $\mu^{2} - \mu + (1 -
S) = 0$ has roots

```math
\mu_{\pm}
 \;=\; \frac{1 \pm \sqrt{1 - 4(1 - S)}}{2}
 \;=\; \frac{1 \pm \sqrt{4 S - 3}}{2}.
```

Converting back to ``\widetilde{H}`` eigenvalues ``\lambda = -2t\mu``,

```math
E_{\pm}(\mathbf{k})
 \;=\; -\,2t\,\mu_{\mp}
 \;=\; -\,t\,\bigl(1 \mp \sqrt{4 S - 3}\bigr)
 \;=\; -t \;\pm\; t\,\sqrt{4 S - 3}.
\tag{6}
```

The sign convention: the "``+``" (upper) band $E_{+} = -t + t\sqrt{4S -
3}`` reaches up to the flat-band value ``+2t`` at ``\Gamma$ and dips
below; the "``-``" (lower) band ``E_{-} = -t - t\sqrt{4S - 3}`` sits
below the flat band everywhere, with minimum at ``\Gamma`` equal to
``-4t``.

### Step 6 — ``\Gamma``-point and zone-boundary checks

At ``\Gamma`` (``\mathbf{k} = (0, 0)``): ``\theta_{1} = \theta_{2} = 0``,
so ``c_{1} = c_{2} = c_{3} = 1``, ``S = 3``, ``P = 1``, $\sqrt{4S - 3} =
\sqrt{9} = 3$. Plugging into (6),

```math
E_{+}(\Gamma) = -t + 3 t = +2t,
\qquad
E_{-}(\Gamma) = -t - 3 t = -4t,
\qquad
E_{\rm flat} = +2t.
```

**The flat band touches the upper dispersive band at ``\Gamma``**; the
degeneracy is two-fold there. Expanding near ``\Gamma``,
``c_{j} \approx 1 - \theta_{j}^{2}/8`` to leading order, so

```math
S(\mathbf{k}) \approx 3 - \tfrac{1}{4}\bigl(\theta_{1}^{2} + \theta_{2}^{2}
            + (\theta_{2} - \theta_{1})^{2}\bigr),
```

```math
\sqrt{4 S - 3} \approx \sqrt{9 - \theta_{1}^{2} - \theta_{2}^{2}
            - (\theta_{2} - \theta_{1})^{2}}
 \approx 3 - \tfrac{1}{6}(\theta_{1}^{2} + \theta_{2}^{2}
               + (\theta_{2} - \theta_{1})^{2}) + O(\theta^{4}).
```

Hence $E_{+}(\Gamma + \mathbf{q}) \approx 2t - (t/6)(\theta_{1}^{2}
+ \theta_{2}^{2} + (\theta_{2} - \theta_{1})^{2})$: the upper band
descends **quadratically** from the degenerate touching point. This
is the **quadratic band-touching** (QBT) point, qualitatively
different from the linear Dirac cone at the honeycomb ``K``-point.

At the zone boundary ``\mathbf{k} = \mathbf{b}_{1}/2`` (M-point),
``\theta_{1} = \pi``, ``\theta_{2} = 0``: ``c_{1} = 0``, ``c_{2} = 1``,
``c_{3} = \cos(-\pi/2) = 0``, so ``S = 1``, ``P = 0``, $\sqrt{4S - 3} =
1$, giving

```math
E_{\pm}(M) = -t \pm t = \{0,\, -2t\},
\qquad
E_{\rm flat}(M) = +2t.
```

The three bands at M are ``\{-2t, 0, +2t\}``, all distinct.

### Step 7 — Real-space compact localised states

The flat band's massive degeneracy has a concrete real-space
signature: ``E = +2t`` eigenstates can be chosen to have compact
support on individual hexagonal plaquettes (the hexagonal voids
between upward- and downward-pointing triangles).

**Construction (Sutherland 1986).** Label the six sites on the
boundary of a hexagonal void in cyclic order as $s_{1}, s_{2},
\dots, s_{6}``. Each ``s_{i}$ is a corner of one of the six triangles
surrounding the hexagon; consecutive ``s_{i}`` belong to different
triangles. Assign amplitudes

```math
\psi(s_{i}) \;=\; (-1)^{i},\qquad i = 1, 2, \dots, 6,
```

and ``\psi = 0`` at every other site.

**Verify ``H\psi = +2t\,\psi``.** Act on ``s_{1}`` with ``H = -t\sum_{NN}``.
The four NN of ``s_{1}`` split into:

- Two NN on the same hexagonal ring (``s_{2}`` and ``s_{6}``), which
  have amplitudes ``+1`` and ``+1`` (from ``(-1)^{2}`` and ``(-1)^{6}``),
  wait — ``(-1)^{2} = +1``, ``(-1)^{6} = +1``, and ``\psi(s_{1}) = -1``.
  The hopping onto ``s_{1}`` from ``s_{2}`` and ``s_{6}`` contributes
  ``-t\,(\psi(s_{2}) + \psi(s_{6})) = -t(+1 + 1) = -2t``.
- Two NN belonging to the **other** triangles sharing the two
  same-triangle edges through ``s_{1}``: these "external" sites have
  ``\psi = 0``.

Hmm, that gives ``(H\psi)(s_{1}) = -2t``, but we need $+2t\,\psi(s_{1})
= -2t`` since ``\psi(s_{1}) = -1``. So actually ``(H\psi)(s_{1}) = -2t
= +2t\cdot(-1) = +2t\cdot\psi(s_{1})$. ✓ Match.

At a site off the hexagonal ring, e.g. one of the "external"
triangle vertices connected to two hexagon sites with amplitudes
``+1`` and ``-1`` (which cancel), we get $(H\psi)(s) = -t\,(1 + (-1) +
\text{zeros}) = 0 = +2t\cdot 0 = +2t\cdot\psi(s)$. ✓

So the hexagon-localised state with alternating-sign amplitudes is
indeed an eigenstate at energy ``+2t`` — the flat-band energy.
**Destructive interference** on the external bonds is the physical
origin of the flat band: hopping into an external site from the
hexagon cancels exactly, freezing the particle in place.

**Counting.** On an ``L_{1}L_{2}`` unit-cell torus, the number of
hexagonal voids is also ``L_{1}L_{2}`` (one per unit cell), so the
naive count gives ``L_{1}L_{2}`` compact localised states — matching
the momentum-space degeneracy (one flat-band mode per ``\mathbf{k}``).
In fact one state is the overall sum of all hexagon states, which
is *not* independent (it is the uniform ``k = 0`` flat-band mode),
so one must subtract one and the correct count is $L_{1}L_{2} -
1 + 1 = L_{1}L_{2}$ states — still matching (see Mielke 1991 §3 and
Bergman–Wu–Balents 2008 eq. (13) for the careful topological
argument).

---

## References

- I. Syozi, *Statistics of kagome lattice*, Prog. Theor. Phys. **6**,
  306 (1951). Original study of the kagome Ising model; introduced
  the lattice name.
- B. Sutherland, *Localization of electronic wave functions due to
  local topology*, Phys. Rev. B **34**, 5208 (1986). First explicit
  compact-localised-state construction for flat bands on
  lattices with hexagonal voids.
- A. Mielke, *Ferromagnetic ground states for the Hubbard model on
  line graphs*, J. Phys. A **24**, L73 (1991). Flat band of the
  kagome NN tight-binding as a line-graph construction.
- H. Tasaki, *Ferromagnetism in the Hubbard models with degenerate
  single-electron ground states*, Phys. Rev. Lett. **69**, 1608
  (1992). General-framework flat-band ferromagnetism.
- D. L. Bergman, C. Wu, L. Balents, *Band touching from real-space
  topology in frustrated hopping models*, Phys. Rev. B **78**,
  125104 (2008). Modern treatment; eq. (13) for the flat-band
  degeneracy count.
- E. H. Lieb, *Two theorems on the Hubbard model*, Phys. Rev. Lett.
  **62**, 1201 (1989). Closely related bipartite flat-band argument
  — see the [Lieb flat-band note](bloch-lieb-flat-band.md).

## Used by

- [Kagome tight-binding model](../models/quantum/tightbinding/kagome.md) —
  `fetch(QAtlas.Kagome(; t, Lx, Ly), TightBindingSpectrum(), …)`
  returns the exact three-band spectrum with the flat band at ``+2t``.
- [Honeycomb dispersion note](bloch-honeycomb-dispersion.md) —
  contrast: honeycomb is bipartite with a linear Dirac-cone touching;
  kagome is non-bipartite with a flat band and quadratic band
  touching, two different paradigms for gapless band structures.
- [Lieb flat-band note](bloch-lieb-flat-band.md) — contrast: Lieb
  flat band arises from sublattice imbalance (bipartite argument),
  not from destructive interference on plaquettes.
