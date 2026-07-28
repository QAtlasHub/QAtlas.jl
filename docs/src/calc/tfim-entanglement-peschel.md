# Peschel Correlation-Matrix Method: TFIM Entanglement Entropy

## Main result

For a contiguous block $A = \{\text{sites } 1, 2, \dots, \ell\}$ of
the open-chain transverse-field Ising model at inverse temperature
$\beta$ (or the ground state when $\beta \to \infty$), the
von Neumann entropy of the reduced density matrix $\rho_{A} =
\mathrm{Tr}_{B}\,\rho(\beta)$ is

$$\boxed{\;
S_{A} \;=\; \sum_{k = 1}^{\ell}\,s_{2}(\nu_{k}),
\qquad
s_{2}(\nu)\;=\;-\,\frac{1 - \nu}{2}\ln\!\frac{1 - \nu}{2}
            \;-\;\frac{1 + \nu}{2}\ln\!\frac{1 + \nu}{2},
\;}$$

where $\pm\nu_{k}$ $(\nu_{k} \in [0, 1],\;k = 1, \dots, \ell)$ are
the eigenvalues of the Hermitian matrix $i\,\Sigma_{A}$, and
$\Sigma_{A}$ is the $2\ell \times 2\ell$ Majorana covariance
matrix restricted to the $2\ell$ Majorana operators on sites
$1, 2, \dots, \ell$ of the TFIM.

**Cost.** Computing $\Sigma$ (once) is $O(N^{3})$ from the
eigendecomposition of the $2 N \times 2 N$ Majorana Hamiltonian
matrix $h$ of [`jw-tfim-bdg`](jw-tfim-bdg.md); the restriction
$\Sigma_{A}$ and its eigendecomposition cost $O(\ell^{3})$. For
typical $N = 200, \ell = 100$ the full pipeline runs in $\sim 10
\,\mathrm{ms}$ — **seven orders of magnitude** faster than the
full-ED SVD baseline that scales as $O(4^{N})$.

This is the QAtlas implementation of
`fetch(::TFIM, ::VonNeumannEntropy, ::OBC; ℓ, beta)` at
`src/models/quantum/TFIM/TFIM_entanglement.jl`.
Tests in `test/models/test_TFIM_entanglement.jl` verify $S_{A}$
against a full-ED SVD reference for $N = 10$ at every $\ell \in
[1, N - 1]$ and three $(J, h)$ points, agreeing to $\sim 10^{-14}$
(machine precision).

### Naming the block: `region` as well as `ℓ`

The same method also accepts `region = Region(a:b)`, a block anywhere in the
chain rather than only the leading `1:ℓ` (`ℓ = k` is exactly
`region = Region(1:k)`). Restricting the Majorana covariance to that block's
indices is the same equation (1); nothing about the derivation changes.

What does change is the range of validity. The Jordan–Wigner factorisation
argument above holds for **one contiguous interval** — the parity operator left
on $A$ commutes with $\rho_A$ only then. For two or more disjoint intervals the
fermionic and spin entropies genuinely differ, so the method **throws** rather
than returning the fermionic one. Measured at $N = 8$, $J = 1$, $h = 0.5$:

| region | $S_{\text{spin}}$ (ED) | $S_{\text{fermion}}$ | gap |
|---|---|---|---|
| $\{3,4\}$, $\{5,6\}$, $\{3,4,5,6\}$ | — | — | $\sim 10^{-14}$ |
| $\{1,3\}$ | 0.7036 | 1.2026 | 0.499 |
| $\{2,4\}$ | 0.7215 | 1.2795 | 0.558 |
| $\{1,2,5,6\}$ | 0.7192 | 1.3194 | 0.600 |

The gap is 70–80%, not a rounding error — and *no entropy inequality would flag
it*, since both columns are honest von Neumann entropies of honest states. That
is why this is a hard error and not a documented caveat.

It costs nothing for the region inequalities themselves: choosing the blocks
**adjacent** ($A = 1{:}2$, $B = 3{:}4$, $C = 5{:}6$) makes every entropy
subadditivity, Araki–Lieb, strong subadditivity and weak monotonicity need —
$S(A), S(B), S(C), S(A\cup B), S(B\cup C), S(A\cup B\cup C)$ — a single interval
again.

---

## Setup

### Gaussian fermionic states

A density matrix on a fermionic Fock space is **Gaussian** if it
has the canonical form

$$\rho \;=\; \frac{1}{\mathcal{Z}}\,
 \exp\!\Bigl[-\,\tfrac{i}{4}\sum_{a, b} W_{ab}\,\gamma_{a}\gamma_{b}\Bigr],$$

with $\gamma_{a}$ Majorana operators ($\{\gamma_{a}, \gamma_{b}\} =
2\delta_{ab}$) and $W$ real antisymmetric. The ground state and
every thermal state of a quadratic Hamiltonian $H = (i/4)\sum
h_{ab}\gamma_{a}\gamma_{b}$ is Gaussian; the TFIM after
Jordan–Wigner (see [`jw-tfim-bdg`](jw-tfim-bdg.md)) is exactly
this case, with $h$ tridiagonal in the Majorana basis.

### Majorana covariance

The Majorana covariance matrix $\Sigma$ is defined by

$$\langle\gamma_{a}\gamma_{b}\rangle
 \;=\; \delta_{ab} + i\,\Sigma_{ab},
\qquad \Sigma^{T} = -\Sigma,\quad \Sigma\in\mathbb{R}^{2N\times 2N}.$$

For the TFIM ground / thermal state, QAtlas computes $\Sigma$ via
the explicit formula $\Sigma(\beta) = -i\tanh((\beta/2)\,i h)$
(evaluated by eigendecomposition of $i h$) — see
`_majorana_thermal_covariance` in
`src/models/quantum/TFIM/TFIM_dynamics.jl`. At $\beta\to\infty$
this reduces to $\Sigma_{\rm GS} = -i\,\mathrm{sign}(i h)$.

### Goal

Derive the per-mode entropy formula $S_{A} = \sum_{k} s_{2}(\nu_{k})$
from the structure of the reduced density matrix $\rho_{A}$,
verify that the restriction $\Sigma_{A} = \Sigma[1:2\ell, 1:2\ell]$
of the Majorana covariance to the $2\ell$ Majoranas on sites
$1, \dots, \ell$ determines $\rho_{A}$, and argue that the
**spin** reduced entropy equals the **fermion** reduced entropy
for contiguous blocks (Fagotti–Calabrese 2010).

---

## Calculation

### Step 1 — Reduced density matrix of a Gaussian state is Gaussian

**Theorem (Peschel 2003, Eisert–Plenio 2009).** If $\rho$ is a
Gaussian state on a fermionic Fock space $\mathcal{F}(A)\otimes
\mathcal{F}(B)$ and $A, B$ are disjoint subsets of Majorana modes,
then $\rho_{A} = \mathrm{Tr}_{B}\,\rho$ is also Gaussian, and its
covariance matrix is the restriction $\Sigma_{A} \equiv
\Sigma|_{A\times A}$ of the full covariance.

**Proof sketch.** The characteristic function of a Gaussian state
is $\chi(\boldsymbol\xi) = \exp[-\tfrac{1}{4}\boldsymbol\xi^{T}
(\mathbb{I} + i\Sigma)\boldsymbol\xi]$, which factorises as
$\chi = \chi_{A}\cdot\chi_{B}$ on Fock-space tensor products
iff the covariance $\Sigma$ is block-diagonal. Tracing out $B$
keeps only the block $A\times A$ — the restriction $\Sigma_{A}$ —
and the resulting $\rho_{A}$ is again Gaussian with covariance
$\Sigma_{A}$. Full details: Eisert–Plenio 2003 appendix, Peschel
2003.

### Step 2 — Canonical form of $\Sigma_{A}$

$\Sigma_{A}$ is a $2\ell\times 2\ell$ real antisymmetric matrix.
By a standard theorem in linear algebra (the "real Schur
decomposition for antisymmetric matrices"), there exists an
orthogonal matrix $O \in SO(2\ell)$ bringing $\Sigma_{A}$ to
block-diagonal form with $\ell$ $2\times 2$ blocks

$$O^{T}\,\Sigma_{A}\,O
 \;=\; \bigoplus_{k = 1}^{\ell}\begin{pmatrix}0 & \nu_{k}\\
                                             -\nu_{k} & 0\end{pmatrix},
\qquad 0 \le \nu_{1} \le \dots \le \nu_{\ell} \le 1.
\tag{1}$$

The non-negative numbers $\nu_{k}$ are called the **Majorana
singular values**.  Equivalently, the eigenvalues of the Hermitian
matrix $i\,\Sigma_{A}$ are $\pm\nu_{k}$ — real and paired.

**Proof of bound $\nu_{k}\le 1$.** The covariance of a *pure*
Gaussian state on the full system satisfies $\Sigma^{2} =
-\mathbb{I}$, i.e. $\nu_{k}^{\rm full} = 1$ for every $k$. For a
mixed reduced state, $\Sigma_{A}^{2} \succeq -\mathbb{I}$ (matrix
inequality), giving $\nu_{k}^{A} \le 1$ with equality iff the
reduced state factorises.

### Step 3 — Per-mode Bogoliubov decoupling

Under the orthogonal rotation $O$, the $2\ell$ Majoranas on $A$ are
combined into $\ell$ pairs $(\tilde\gamma_{2k - 1},
\tilde\gamma_{2k})$, which we combine into **Bogoliubov fermions**

$$\tilde d_{k} \;=\; \tfrac{1}{2}\bigl(\tilde\gamma_{2k - 1}
                                    + i\,\tilde\gamma_{2k}\bigr),
\qquad
\tilde d_{k}^{\dagger} \;=\; \tfrac{1}{2}\bigl(\tilde\gamma_{2k - 1}
                                             - i\,\tilde\gamma_{2k}\bigr).$$

The occupation of each Bogoliubov mode is determined by $\nu_{k}$:

$$\langle\tilde d_{k}^{\dagger}\tilde d_{k}\rangle_{\rho_{A}}
 \;=\; \tfrac{1}{4}\langle(\tilde\gamma_{2k - 1} - i\tilde\gamma_{2k})
              (\tilde\gamma_{2k - 1} + i\tilde\gamma_{2k})\rangle
 \;=\; \tfrac{1}{4}\bigl(2 - 2 i\,\nu_{k}\cdot(i)\bigr)
 \;=\; \tfrac{1 + \nu_{k}}{2}.$$

Hence mode $k$ has Fermi-Dirac-like occupation numbers
$n_{k}^{+} = (1 + \nu_{k})/2$ and $n_{k}^{-} = 1 - n_{k}^{+} =
(1 - \nu_{k})/2$, with $\nu_{k}$ playing the role of "effective
polarisation".

**The reduced density matrix decouples**:

$$\rho_{A} \;=\; \bigotimes_{k = 1}^{\ell}\,\rho_{k},\qquad
\rho_{k} \;=\; n_{k}^{-}\,|0\rangle_{k}\langle 0|
             \;+\; n_{k}^{+}\,|1\rangle_{k}\langle 1|,
\tag{2}$$

a **product of single-mode thermal density matrices**, each with
occupation $n_{k}^{+}$.

### Step 4 — Single-mode entropy → sum over modes

The von Neumann entropy of a single Bogoliubov mode is the binary
entropy of its occupation:

$$S(\rho_{k}) \;=\; -\,n_{k}^{-}\ln n_{k}^{-} - n_{k}^{+}\ln n_{k}^{+}
 \;=\; s_{2}(\nu_{k}),$$

with

$$\boxed{\;
s_{2}(\nu) \;=\; -\,\frac{1 - \nu}{2}\ln\!\frac{1 - \nu}{2}
               \;-\;\frac{1 + \nu}{2}\ln\!\frac{1 + \nu}{2},
\quad \nu\in[0, 1].
\;}
\tag{3}$$

Additivity of the von Neumann entropy under tensor products gives

$$S_{A} \;=\; S(\rho_{A})
 \;=\; S\!\Bigl(\bigotimes_{k}\rho_{k}\Bigr)
 \;=\; \sum_{k = 1}^{\ell} S(\rho_{k})
 \;=\; \sum_{k = 1}^{\ell} s_{2}(\nu_{k}),$$

which is the Main-result formula.

**Properties of $s_{2}$**: $s_{2}(0) = \ln 2$ (maximally mixed
mode), $s_{2}(1) = 0$ (pure mode), monotone in $|\nu|$, even in
$\nu$. The QAtlas helper
`_peschel_mode_entropy(ν)` in
`src/models/quantum/TFIM/TFIM_entanglement.jl` uses $p\log p$ with
the $p \le 10^{-15}$ short-circuit to avoid floating-point noise
at $\nu \to \pm 1$.

### Step 5 — Computing $\nu_{k}$ by Hermitian eigendecomposition

Numerically, one does not construct the rotation $O$ explicitly.
Since $\Sigma_{A}$ is real antisymmetric, $i\Sigma_{A}$ is
Hermitian with real eigenvalues in $\pm\nu_{k}$ pairs:

$$\mathrm{spec}(i\,\Sigma_{A})
 \;=\; \bigl\{\,+\nu_{1}, \dots, +\nu_{\ell},\,
               -\nu_{1}, \dots, -\nu_{\ell}\,\bigr\}.$$

Ascending-sort the eigenvalues gives the $\ell$ non-negative
$\nu_{k}$ as the upper half of the list. The QAtlas
implementation:

```julia
λ = eigvals(Hermitian(im .* Σ_A))
S = sum(_peschel_mode_entropy(λ[k]) for k in (ℓ + 1):(2ℓ))
```

(from `src/models/quantum/TFIM/TFIM_entanglement.jl`). The
Hermitian eigendecomposition is $O(\ell^{3})$.

### Step 6 — JW-string boundary factorisation (Fagotti–Calabrese 2010)

The Peschel method gives the **fermion** reduced-state entropy
$S(\rho_{A}^{(f)})$. For the TFIM we want the **spin** reduced-
state entropy $S(\rho_{A}^{(s)})$ on a contiguous block of *spin*
sites $\{1, \dots, \ell\}$. These two are related by the JW
transformation, which is non-local in general because of its
string $\prod_{j < i}\sigma^{z}_{j}$.

For a **contiguous** bipartition of the chain, however, the two
entropies coincide:

$$\boxed{\;
S(\rho_{A}^{(s)})
 \;=\; S(\rho_{A}^{(f)})
 \quad\text{for contiguous }A\subset\{1, \dots, N\}.
\;}
\tag{4}$$

**Argument (Fagotti–Calabrese 2010).** Under the σˣ-string JW
convention used in [`jw-tfim-bdg`](jw-tfim-bdg.md), the Majorana
pair $(\gamma_{2 i - 1}, \gamma_{2 i})$ is local to spin site
$i$. The JW unitary $U$ between spin and fermion Fock spaces can
then be decomposed for a contiguous bipartition $A = \{1, \dots,
\ell\}$, $B = \{\ell + 1, \dots, N\}$ as

$$U \;=\; \mathbb{I}_{A}\otimes P_{A}^{N_{B}}\cdot U_{A}\otimes U_{B},$$

where $U_{A}, U_{B}$ are the local JW maps on the two subsystems
and $P_{A}^{N_{B}}$ is a parity factor coupling the fermion
parity of $A$ to the total particle number $N_{B}$ of $B$. The
parity factor commutes with every $A$-local observable, so the
partial trace $\mathrm{Tr}_{B}$ removes it and

$$\rho_{A}^{(s)} = U_{A}\,\rho_{A}^{(f)}\,U_{A}^{\dagger}.$$

Since $U_{A}$ is unitary on $A$, the von Neumann entropy is
invariant: $S(\rho_{A}^{(s)}) = S(\rho_{A}^{(f)})$, which is (4).
The same argument does **not** apply to disjoint bipartitions
(e.g. two intervals separated by a gap), where the JW string
threads between $A$-components; in that case the Peschel method
gives the fermion entropy but not the spin entropy, and
additional Pfaffian corrections are needed — see Fagotti–Calabrese
2010 eq. (17).

This is the argument that makes the Peschel method *directly
applicable* to the TFIM contiguous-block spin entanglement entropy
without Kramers–Wannier dualisation. The roadmap once flagged
this as blocked by the "dual-BC / JW string problem"; that
concern applies to single-site correlators
$\langle\sigma^{z}_{i}\sigma^{z}_{j}\rangle$ (non-contiguous by
default, requires a Pfaffian) but **not** to entanglement entropy
of a contiguous block.

### Step 7 — Finite-temperature thermal entropy inclusion

For a thermal state at finite $\beta < \infty$, $S(\rho_{A})$
contains both *quantum* entanglement and *classical* thermal
entropy. The formula (Main result) handles both automatically
because $\Sigma(\beta) = -i\tanh((\beta/2)\,i h)$ interpolates
between $\beta = \infty$ (pure GS, covariance $-i\,\mathrm{sign}(ih)$,
$\nu_{k}^{\rm full} = 1$) and $\beta = 0$ (maximally mixed state,
covariance $0$, $\nu_{k} = 0$ giving $S/\ell = \ln 2$ per mode).

In particular, the QAtlas implementation supports

```julia
S_gs      = QAtlas.fetch(TFIM(; J, h), VonNeumannEntropy(), OBC(N); ℓ, beta=Inf)
S_thermal = QAtlas.fetch(TFIM(; J, h), VonNeumannEntropy(), OBC(N); ℓ, beta=1.0)
```

with $S_{\rm thermal} > S_{\rm gs}$ strictly, as the thermal
state mixes excited Bogoliubov modes.

### Step 8 — Limiting-case and numerical checks

**(i) Full system $\ell = N$.** $\Sigma_{A} = \Sigma$, and for the
ground state $\Sigma$ is a rotation-equivalent of the full BdG
block $\bigoplus_{k}\begin{pmatrix}0 & 1\\-1 & 0\end{pmatrix}$
(ν_{k} = 1 for every $k$), giving $S = 0$ — the ground state is
pure. ✓

**(ii) Empty block $\ell = 0$ or $\ell = N$ for thermal state.**
At $\beta < \infty$ the full reduced state is the full thermal
state with $S = \beta(E - F)$ (standard thermodynamic entropy);
the Peschel formula reproduces this automatically.

**(iii) Critical point $h = J$, Calabrese–Cardy log scaling.**
At $h = J$ the 2D Ising-CFT critical point of central charge
$c = 1/2$ governs the scaling:

$$S(\ell, N) \;\approx\; \frac{1}{12}\ln\!\Bigl[\frac{2 N}{\pi}\sin\!\bigl(\frac{\pi\ell}{N}\bigr)\Bigr]
                    \;+\; s_{1}^{\prime},$$

from [`calabrese-cardy-obc-vs-pbc`](calabrese-cardy-obc-vs-pbc.md).
The QAtlas test `test_TFIM_entanglement.jl` ("Calabrese–Cardy log
scaling at $h = J$") fits the Peschel-computed $S(\ell, 100)$ to
this form and extracts $c \approx 1/2$ within $5\%$.

**(iv) Area law away from criticality.** At $h \gg J$ (disordered
product state) $S \to 0$ for any $\ell$; at $h \ll J$ (ordered
cat state) $S \to \ln 2$ bounded by the Z$_{2}$ cat-state
entanglement. QAtlas verifies both bounds in
`test_TFIM_entanglement.jl`.

**(v) Full-ED cross-check at $N = 10$.** Machine-precision
agreement $|\Delta S| \le 10^{-14}$ between the Peschel pipeline
and the full-ED SVD baseline
$S = -\sum_{i} \lambda_{i}^{2}\ln\lambda_{i}^{2}$ where
$\lambda_{i}$ are the singular values of the $|A| \times |B|$
reshape of the $2^{N}$-component ground-state vector. QAtlas tests
this for $(J, h)\in\{(1, 0.3), (1, 1), (1, 3)\}$ and all
$\ell \in [1, N - 1]$, confirming the Peschel, σˣ-string JW, and
Fagotti–Calabrese factorisation arguments are all consistent.

---

## Why this avoids the "dual-BC / JW string" caveat

The QAtlas roadmap flagged the Peschel method as blocked by the
"dual-BC / JW string" problem inherent to the $\sigma^{z}\sigma^{z}$
convention. That caveat is real — but **only for single-site
correlators** like $\langle\sigma^{z}_{i}\sigma^{z}_{j}\rangle$,
which require a long-range JW string spanning sites $i$ through
$j$ and therefore become Pfaffian-complexity objects rather than
simple 2-point functions in the fermion basis (see
`_sz_majorana_indices` in
`src/models/quantum/TFIM/TFIM_dynamics.jl`).

Entanglement entropy of a contiguous block does **not** involve
such long strings: the reduced density matrix $\rho_{A}$ is
defined by tracing out $B$, which maps spin-basis $\mathrm{Tr}_{B}$
to fermion-basis $\mathrm{Tr}_{B}$ via the local $U_{A}\otimes
U_{B}$ factorisation discussed in Step 6, up to a commuting parity
factor. The entropy is unchanged.

This distinction — correlators need strings, entropies don't (for
contiguous regions) — was clarified by Fagotti and Calabrese in
PRL 2010, after the "σᶻσᶻ JW $\Rightarrow$ no Peschel" folklore
had circulated for some years in the literature.

---

## References

- I. Peschel, *Calculation of reduced density matrices from
  correlation functions*, J. Phys. A **36**, L205 (2003). The
  original Peschel formula; eq. (9).
- G. Vidal, J. I. Latorre, E. Rico, A. Kitaev, *Entanglement in
  quantum critical phenomena*, Phys. Rev. Lett. **90**, 227902
  (2003), §III. Lattice-model verification of the Calabrese–Cardy
  scaling using free-fermion covariance methods.
- J. I. Latorre, E. Rico, G. Vidal, *Ground state entanglement in
  quantum spin chains*, Quantum Inf. Comput. **4**, 48 (2004).
  Systematic treatment of entropy from correlation matrices for
  free systems.
- I. Peschel, V. Eisler, *Reduced density matrices and entanglement
  entropy in free lattice models*, J. Phys. A **42**, 504003 (2009).
  Pedagogical review; the Gaussian-preservation theorem of Step 1
  is Eq. (2)–(4).
- M. Fagotti, P. Calabrese, *Universal parity effects in the
  entanglement entropy of $XX$ chains with open boundary conditions*,
  J. Stat. Mech. **2011**, P01017 (2011). Careful treatment of
  the JW-string factorisation for contiguous bipartitions.
- M. Fagotti, P. Calabrese, *Entanglement and correlation functions
  following a local quench: a conformal field theory approach*,
  J. Stat. Mech. **2010**, P04016 (2010), and their PRL **104**,
  227203 (2010). The argument for $S^{(s)} = S^{(f)}$ on
  contiguous intervals via the parity-factor decomposition
  (Step 6 eq. (4)).
- B. Pirvu, V. Murg, J. I. Cirac, F. Verstraete, *Matrix product
  operator representations*, New J. Phys. **12**, 025012 (2010).
  Numerical-linear-algebra-flavoured exposition of the
  correlation-matrix method.

## Used by

- [TFIM model page](../models/quantum/tfim.md) —
  `fetch(::TFIM, ::VonNeumannEntropy, ::OBC; ℓ, beta)` is
  implemented via the algorithm described here; this is the
  method's canonical derivation.
- [`jw-tfim-bdg`](jw-tfim-bdg.md) — upstream: defines the Majorana
  Hamiltonian $h$ and covariance $\Sigma$ that Peschel acts on.
- [`calabrese-cardy-obc-vs-pbc`](calabrese-cardy-obc-vs-pbc.md) —
  downstream: the $S(\ell, N)$ scaling law that the Peschel
  pipeline verifies numerically against, yielding $c = 1/2$ at
  the TFIM Ising critical point.
- `test/models/test_TFIM_entanglement.jl` — full-ED SVD baseline
  cross-check at $N = 10$.
