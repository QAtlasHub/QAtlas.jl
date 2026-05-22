# Heisenberg Dimer: Singlet–Triplet Splitting

## Main result

For two spin-``\tfrac{1}{2}`` particles coupled by the isotropic
Heisenberg exchange

```math
H \;=\; J\,\mathbf{S}_{1}\cdot\mathbf{S}_{2},
```

the Hilbert space ``\mathbb{C}^{2}\otimes\mathbb{C}^{2}`` (dimension
``4``) decomposes into a singlet (``S_{\rm tot} = 0``) and a triplet
(``S_{\rm tot} = 1``), with spectrum and gap

```math
\boxed{\;
\mathrm{Spectrum}(H) \;=\; \Bigl\{\,-\tfrac{3 J}{4},\;\tfrac{J}{4},\;
                                \tfrac{J}{4},\;\tfrac{J}{4}\,\Bigr\},
\qquad
\Delta \;=\; E_{t} - E_{s} \;=\; J.
\;}
```

For ``J > 0`` (antiferromagnetic) the singlet is the ground state;
for ``J < 0`` (ferromagnetic) the three degenerate triplet states
are the ground manifold.

The singlet and triplet eigenvectors are

```math
\boxed{\;
|s\rangle = \tfrac{1}{\sqrt{2}}\bigl(|{\uparrow\downarrow}\rangle - |{\downarrow\uparrow}\rangle\bigr),
\qquad
\{|t_{+1}\rangle,\,|t_{0}\rangle,\,|t_{-1}\rangle\}
 = \Bigl\{|{\uparrow\uparrow}\rangle,\;\tfrac{1}{\sqrt{2}}\bigl(|{\uparrow\downarrow}\rangle + |{\downarrow\uparrow}\rangle\bigr),\;|{\downarrow\downarrow}\rangle\Bigr\}.
\;}
```

This is the ``N = 2`` building block of every result for spin-1/2
Heisenberg chains and ladders; its exact spectrum is the reference
point against which every lattice-ED and QAtlas cross-check is
validated at small system sizes.

---

## Setup

### Operators and conventions

Spin-1/2 operators $\mathbf{S}_{i} = \tfrac{1}{2}
\boldsymbol{\sigma}_{i}`` act on the ``i``-th factor of ``\mathcal{H} =
\mathbb{C}^{2}\otimes\mathbb{C}^{2}$. The computational basis is

```math
\bigl\{\,|{\uparrow\uparrow}\rangle,\,
         |{\uparrow\downarrow}\rangle,\,
         |{\downarrow\uparrow}\rangle,\,
         |{\downarrow\downarrow}\rangle\,\bigr\}.
```

Raising / lowering operators
``S^{\pm}_{i} = S^{x}_{i} \pm i S^{y}_{i}`` satisfy
``S^{+}_{i}|{\downarrow}\rangle_{i} = |{\uparrow}\rangle_{i}``,
``S^{-}_{i}|{\uparrow}\rangle_{i} = |{\downarrow}\rangle_{i}``, all
other actions zero. The Casimir

```math
\mathbf{S}_{i}^{2} \;=\; \tfrac{3}{4}\,\mathbb{I}_{i},
\qquad s = \tfrac{1}{2},\ \ s(s + 1) = \tfrac{3}{4}.
```

### Total-spin decomposition

Angular-momentum addition ``1/2 \otimes 1/2 = 0 \oplus 1`` splits the
four basis states into a singlet and a triplet. We will label
eigenstates by the two commuting good quantum numbers
``(S_{\rm tot}, S^{z}_{\rm tot})``.

### Goal

Derive the four eigenvalues of ``H`` from the ``\mathbf{S}_{\rm tot}^{2}``
Casimir, identify the corresponding eigenvectors explicitly, and
verify them by direct action of ``H`` and ``\mathbf{S}^{2}_{\rm tot}``
on the basis.

---

## Calculation

### Step 1 — Rewrite ``H`` using ``\mathbf{S}_{\rm tot}^{2}``

Define ``\mathbf{S}_{\rm tot} = \mathbf{S}_{1} + \mathbf{S}_{2}`` and
square:

```math
\mathbf{S}_{\rm tot}^{2}
 \;=\; (\mathbf{S}_{1} + \mathbf{S}_{2})^{2}
 \;=\; \mathbf{S}_{1}^{2} + \mathbf{S}_{2}^{2} + 2\,\mathbf{S}_{1}\cdot\mathbf{S}_{2}.
```

Solve for ``\mathbf{S}_{1}\cdot\mathbf{S}_{2}``:

```math
\mathbf{S}_{1}\cdot\mathbf{S}_{2}
 \;=\; \tfrac{1}{2}\bigl[\mathbf{S}_{\rm tot}^{2}
                        - \mathbf{S}_{1}^{2} - \mathbf{S}_{2}^{2}\bigr]
 \;=\; \tfrac{1}{2}\bigl[\,S_{\rm tot}(S_{\rm tot} + 1) - \tfrac{3}{2}\,\bigr].
\tag{1}
```

Multiplying by ``J``, we read off the two eigenvalues of ``H``:

```math
E(S_{\rm tot} = 0)
 = J\cdot\tfrac{1}{2}\bigl[\,0 - \tfrac{3}{2}\,\bigr]
 = -\tfrac{3 J}{4},
```

```math
E(S_{\rm tot} = 1)
 = J\cdot\tfrac{1}{2}\bigl[\,2 - \tfrac{3}{2}\,\bigr]
 = \tfrac{J}{4}.
```

The singlet is one-dimensional (``2S_{\rm tot} + 1 = 1`` state) and
the triplet is three-dimensional (``2S_{\rm tot} + 1 = 3`` states),
matching the Hilbert-space dimension ``1 + 3 = 4``. The **gap** is

```math
\Delta \;=\; E(1) - E(0) \;=\; \tfrac{J}{4} - \bigl(-\tfrac{3 J}{4}\bigr)
 \;=\; J.
```

### Step 2 — Singlet eigenvector by explicit ``\mathbf{S}^{2}_{\rm tot}`` verification

Claim: $|s\rangle = \tfrac{1}{\sqrt{2}}(|{\uparrow\downarrow}\rangle
- |{\downarrow\uparrow}\rangle)`` has ``S_{\rm tot} = 0$ and
``S^{z}_{\rm tot} = 0``.

**Verify ``S^{z}_{\rm tot}|s\rangle = 0``.** Acting on each term:

```math
S^{z}_{\rm tot}|{\uparrow\downarrow}\rangle
 = (\tfrac{1}{2} + (-\tfrac{1}{2}))\,|{\uparrow\downarrow}\rangle = 0,
```

```math
S^{z}_{\rm tot}|{\downarrow\uparrow}\rangle
 = (-\tfrac{1}{2} + \tfrac{1}{2})\,|{\downarrow\uparrow}\rangle = 0.
```

So ``S^{z}_{\rm tot}|s\rangle = 0``. ✓

**Verify ``\mathbf{S}^{2}_{\rm tot}|s\rangle = 0``.** Use
$\mathbf{S}^{2}_{\rm tot} = (S^{z}_{\rm tot})^{2} + \tfrac{1}{2}
(S^{+}_{\rm tot}S^{-}_{\rm tot} + S^{-}_{\rm tot}S^{+}_{\rm tot})$
with ``S^{\pm}_{\rm tot} = S^{\pm}_{1} + S^{\pm}_{2}``:

```math
S^{+}_{\rm tot}|{\uparrow\downarrow}\rangle
 = S^{+}_{1}|{\uparrow\downarrow}\rangle + S^{+}_{2}|{\uparrow\downarrow}\rangle
 = 0 + |{\uparrow\uparrow}\rangle
 = |{\uparrow\uparrow}\rangle,
```

```math
S^{+}_{\rm tot}|{\downarrow\uparrow}\rangle
 = |{\uparrow\uparrow}\rangle + 0 = |{\uparrow\uparrow}\rangle,
```

so $S^{+}_{\rm tot}|s\rangle = \tfrac{1}{\sqrt{2}}
(|{\uparrow\uparrow}\rangle - |{\uparrow\uparrow}\rangle) = 0$.
Similarly ``S^{-}_{\rm tot}|s\rangle = 0``. Combined with
``(S^{z}_{\rm tot})^{2}|s\rangle = 0``, this gives
``\mathbf{S}^{2}_{\rm tot}|s\rangle = 0``, confirming $S_{\rm tot}
(S_{\rm tot} + 1) = 0 \Rightarrow S_{\rm tot} = 0$. ✓

### Step 3 — Triplet eigenvectors by ``S^{\pm}_{\rm tot}`` ladder

The ``|t_{+1}\rangle = |{\uparrow\uparrow}\rangle`` state has
$S^{z}_{\rm tot}|{\uparrow\uparrow}\rangle = +1 \cdot
|{\uparrow\uparrow}\rangle`` (maximal ``S^z$), so it is the
"highest-weight" vector of the triplet. Apply ``S^{-}_{\rm tot}``:

```math
S^{-}_{\rm tot}|{\uparrow\uparrow}\rangle
 = S^{-}_{1}|{\uparrow\uparrow}\rangle + S^{-}_{2}|{\uparrow\uparrow}\rangle
 = |{\downarrow\uparrow}\rangle + |{\uparrow\downarrow}\rangle
 = \sqrt{2}\cdot\tfrac{1}{\sqrt{2}}\bigl(|{\uparrow\downarrow}\rangle
                                        + |{\downarrow\uparrow}\rangle\bigr),
```

using the standard ladder normalisation $S^{-}|S, m\rangle =
\sqrt{S(S + 1) - m(m - 1)}\,|S, m - 1\rangle`` which for ``S = 1,
m = 1`` gives the factor ``\sqrt{2}$. Hence

```math
|t_{0}\rangle \;=\; \tfrac{1}{\sqrt{2}}\bigl(|{\uparrow\downarrow}\rangle
                                              + |{\downarrow\uparrow}\rangle\bigr).
```

Applying ``S^{-}_{\rm tot}`` once more:

```math
S^{-}_{\rm tot}|t_{0}\rangle
 = \tfrac{1}{\sqrt{2}}\bigl(S^{-}_{\rm tot}|{\uparrow\downarrow}\rangle
                          + S^{-}_{\rm tot}|{\downarrow\uparrow}\rangle\bigr)
 = \tfrac{1}{\sqrt{2}}(|{\downarrow\downarrow}\rangle + 0 + 0 + |{\downarrow\downarrow}\rangle)
 = \sqrt{2}\,|{\downarrow\downarrow}\rangle,
```

confirming ``|t_{-1}\rangle = |{\downarrow\downarrow}\rangle``. The
three triplet states form a Clebsch–Gordan-standard spin-1
multiplet.

**Verify triplet is eigenspace with ``\mathbf{S}^{2}_{\rm tot} = 2``.**
For ``|t_{+1}\rangle``, ``S^{z}|t_{+1}\rangle = +1 |t_{+1}\rangle`` and
``S^{+}|t_{+1}\rangle = 0`` (highest weight). Using
``\mathbf{S}^{2} = (S^{z})^{2} + S^{z} + S^{-}S^{+}``:

```math
\mathbf{S}^{2}|t_{+1}\rangle = (1 + 1 + 0)|t_{+1}\rangle = 2\,|t_{+1}\rangle,
```

i.e. ``S_{\rm tot}(S_{\rm tot} + 1) = 2 \Rightarrow S_{\rm tot} = 1``.
The other two triplet states follow by ladder-preservation of the
``\mathbf{S}^{2}`` eigenvalue.

### Step 4 — Direct ``4\times 4`` diagonalisation cross-check

Write ``H = J\,\mathbf{S}_{1}\cdot\mathbf{S}_{2}`` as an explicit $4
\times 4`` matrix in the computational basis ``\{|{\uparrow\uparrow}
\rangle, |{\uparrow\downarrow}\rangle, |{\downarrow\uparrow}
\rangle, |{\downarrow\downarrow}\rangle\}$. Expand

```math
\mathbf{S}_{1}\cdot\mathbf{S}_{2}
 = S^{z}_{1}S^{z}_{2}
 + \tfrac{1}{2}\bigl(S^{+}_{1}S^{-}_{2} + S^{-}_{1}S^{+}_{2}\bigr).
```

Compute matrix elements:

```math
\begin{aligned}
\langle{\uparrow\uparrow}|\mathbf{S}_{1}\cdot\mathbf{S}_{2}|{\uparrow\uparrow}\rangle
&= \tfrac{1}{2}\cdot\tfrac{1}{2} = \tfrac{1}{4},\\
\langle{\uparrow\downarrow}|\mathbf{S}_{1}\cdot\mathbf{S}_{2}|{\uparrow\downarrow}\rangle
&= \tfrac{1}{2}\cdot(-\tfrac{1}{2}) = -\tfrac{1}{4},\\
\langle{\downarrow\uparrow}|\mathbf{S}_{1}\cdot\mathbf{S}_{2}|{\downarrow\uparrow}\rangle
&= -\tfrac{1}{4},\\
\langle{\downarrow\downarrow}|\mathbf{S}_{1}\cdot\mathbf{S}_{2}|{\downarrow\downarrow}\rangle
&= \tfrac{1}{4},\\
\langle{\uparrow\downarrow}|\mathbf{S}_{1}\cdot\mathbf{S}_{2}|{\downarrow\uparrow}\rangle
&= \tfrac{1}{2},\\
\langle{\downarrow\uparrow}|\mathbf{S}_{1}\cdot\mathbf{S}_{2}|{\uparrow\downarrow}\rangle
&= \tfrac{1}{2}.
\end{aligned}
```

The ``4\times 4`` matrix of ``H/J`` (in the order listed above) is

```math
\frac{H}{J} \;=\; \begin{pmatrix}
 1/4 & 0 & 0 & 0 \\
 0 & -1/4 & 1/2 & 0 \\
 0 & 1/2 & -1/4 & 0 \\
 0 & 0 & 0 & 1/4
\end{pmatrix}.
```

The ``(|{\uparrow\uparrow}\rangle, |{\downarrow\downarrow}\rangle)``
states are already eigenvectors with eigenvalue ``1/4``. The
middle ``2\times 2`` block

```math
M \;=\; \begin{pmatrix}-1/4 & 1/2 \\ 1/2 & -1/4\end{pmatrix}
```

has eigenvalues ``-1/4 \pm 1/2 = \{1/4, -3/4\}`` with eigenvectors
``(1, 1)^{T}/\sqrt{2}`` and ``(1, -1)^{T}/\sqrt{2}``. Multiplying by
``J`` recovers

```math
\mathrm{spec}(H) = \{1/4, 1/4, 1/4, -3/4\}\cdot J
 = \{J/4, J/4, J/4, -3J/4\},
```

exactly the Main-result spectrum. The ``(+, +)`` eigenvector is
``|t_{0}\rangle`` at eigenvalue ``J/4``, and the ``(+, -)`` eigenvector is
``|s\rangle`` at eigenvalue ``-3J/4`` — matching the total-spin
identification.

### Step 5 — Limiting-case checks

**(i) AF ground state ``J > 0``.** The singlet at ``E_{s} = -3J/4`` is
the unique ground state. The first-excited triplet at $E_{t} =
J/4`` is three-fold degenerate, with gap ``\Delta = J$.

**(ii) FM ground state ``J < 0``.** The triplet is now the
three-fold-degenerate ground manifold at ``E_{t} = J/4 < 0``, and
the singlet at ``E_{s} = -3J/4 > 0`` is an excited state. The sign
convention is reversed but the spectrum is the same four
eigenvalues up to reordering.

**(iii) Trace.** $\mathrm{Tr}\,H = J\,\mathrm{Tr}\,(\mathbf{S}_{1}
\cdot\mathbf{S}_{2}) = 0`` identically (since ``\mathbf{S}_{i}$ are
traceless). Check: ``-3 J/4 + 3 \cdot J/4 = 0``. ✓

**(iv) Cross-reference with the spin-1/2 Heisenberg chain
thermodynamic limit.** The infinite AF Heisenberg chain has
ground-state energy per bond $e_{0} = J(1/4 - \ln 2) \approx
-0.443\,J$, derived in
[`bethe-ansatz-heisenberg-e0`](bethe-ansatz-heisenberg-e0.md). For
``N = 2`` the per-bond energy is ``E_{s}/1 = -3 J/4 = -0.75\,J`` (only
one bond on a 2-site chain with open boundary condition; for PBC
there are two "bonds" that double-count, giving ``-3 J/2`` total and
``-3 J/4`` per bond after division by 2). The finite-``N = 2`` value
overshoots the thermodynamic limit because the dimer is maximally
"entangled" locally without the spinon-fluctuation dilution that
occurs at larger ``N``.

### Step 6 — QAtlas lattice verification

QAtlas confirms the dimer spectrum against independent lattice ED:

```julia
λ = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=2, J=1.0, bc=:OBC)
# → [-0.75, 0.25, 0.25, 0.25]
```

matching ``\{-3J/4, J/4, J/4, J/4\}`` for ``J = 1``.
`test/verification/test_heisenberg_dimer.jl` cross-checks this against
a Lattice2D-based spin-1/2 ED build — both constructions give
identical eigenvalues to machine precision.

---

## References

- A. Auerbach, *Interacting Electrons and Quantum Magnetism*
  (Springer, 1994), §2. Textbook derivation of (1) and the singlet–
  triplet decomposition.
- J. J. Sakurai, *Modern Quantum Mechanics*, 2nd ed. (Addison-
  Wesley, 2011), Ch. 3. Clebsch–Gordan addition for two spin-1/2
  particles.
- K. Huang, *Statistical Mechanics*, 2nd ed. (Wiley, 1987), §13.2.
  Heisenberg dimer as the simplest interacting magnetic system.

## Used by

- [Heisenberg model page](../models/quantum/heisenberg.md) —
  ``N = 2`` `ExactSpectrum` returns exactly the Main-result set.
- [`bethe-ansatz-heisenberg-e0`](bethe-ansatz-heisenberg-e0.md) —
  the ``N = 2`` dimer result is the smallest-``N`` limiting-case check
  for the thermodynamic-limit Bethe-ansatz energy $e_{0} =
  J(1/4 - \ln 2)``; the dimer overshoots ``e_{0}$ from below at finite
  ``N`` and converges from finite-size corrections.
- [XXZ model page](../models/quantum/xxz.md) — dimer limit at the
  isotropic point ``\Delta = 1``.
