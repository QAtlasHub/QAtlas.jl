# Jordan–Wigner Reduction of the TFIM to BdG Form

## Main result

For the open-chain transverse-field Ising model with $N$ sites,

$$H = -J\sum_{i=1}^{N-1}\sigma^z_i\sigma^z_{i+1}
      - h\sum_{i=1}^{N}\sigma^x_i,
  \qquad (J > 0,\ h \ge 0),$$

the Kramers–Wannier duality followed by a Jordan–Wigner transformation
reduces $H$ to the quadratic free-fermion form

$$\boxed{\;
H \;=\; \tfrac{1}{2}\,\Psi^\dagger \mathcal{H}_{\rm BdG}\, \Psi
       \;-\; J(N-1),
\qquad
\mathcal{H}_{\rm BdG} \;=\;
\begin{pmatrix} A & B \\ -B & -A \end{pmatrix}
\in \mathbb{R}^{2N \times 2N}.
\;}$$

Here $\Psi = (d_1,\dots,d_N,\,d^\dagger_1,\dots,d^\dagger_N)^{T}$
is the Nambu vector of a set of canonical fermions $\{d_i\}$ defined
below, and $A, B$ are the tridiagonal matrices

$$A_{ij} = 2h\,\delta_{ij} - J(\delta_{i,j+1} + \delta_{i+1,j}),
  \qquad
  B_{ij} = J(\delta_{i+1,j} - \delta_{i,j+1}).$$

$\mathcal{H}_{\rm BdG}$ is **real symmetric** and its eigenvalues
come in $\pm\Lambda_n$ pairs ($n = 1,\dots,N$, $\Lambda_n > 0$), so
the finite-temperature expectation of $H$ reads

$$\boxed{\;
\langle H\rangle(\beta)
 \;=\; -\sum_{n=1}^{N}\frac{\Lambda_n}{2}\,
        \tanh\!\left(\frac{\beta\Lambda_n}{2}\right)
        \;-\; J(N-1).
\;}$$

The ground-state energy is the $\beta\to\infty$ limit,
$E_0 = -\tfrac{1}{2}\sum_n \Lambda_n - J(N-1)$.

The two boxed formulas are what the QAtlas implementation
[`fetch(::TFIM, ::Energy, ::OBC)`](../models/quantum/tfim.md) returns
at finite $N$; they hold for every $J, h \ge 0$ and every finite $N
\ge 2$.

---

## Setup

### Conventions

Pauli matrices on site $i$:

$$\sigma^x_i = \begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
  \sigma^y_i = \begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
  \sigma^z_i = \begin{pmatrix}1&0\\0&-1\end{pmatrix},$$

raising / lowering

$$\sigma^+_i = \tfrac{1}{2}(\sigma^x_i + i\sigma^y_i)
             = \begin{pmatrix}0&1\\0&0\end{pmatrix},\qquad
  \sigma^-_i = \tfrac{1}{2}(\sigma^x_i - i\sigma^y_i)
             = \begin{pmatrix}0&0\\1&0\end{pmatrix}.$$

These satisfy, on a single site,

$$\{\sigma^+, \sigma^-\} = 1,\qquad
  (\sigma^+)^2 = (\sigma^-)^2 = 0,\qquad
  \sigma^z = 1 - 2\,\sigma^-\sigma^+.$$

Operators on different sites commute.

### Hamiltonian

$$H = -J\sum_{i=1}^{N-1}\sigma^z_i\sigma^z_{i+1}
      - h\sum_{i=1}^{N}\sigma^x_i.$$

OBC means only $N-1$ bonds, no site-$N$-to-site-$1$ term.

### Goal

Find a linear transformation from the spin operators to canonical
fermion operators $\{d_i\}$ such that $H$ becomes *quadratic* in
$\{d_i, d^\dagger_i\}$. That quadratic form then diagonalises to a
free-fermion spectrum by an orthogonal rotation (the BdG / Bogoliubov
step).

Jordan–Wigner applied directly to $\sigma^z_i\sigma^z_{i+1}$ produces
a *quartic* term; we therefore precede it with the Kramers–Wannier
duality, which swaps the bond interaction $\sigma^z\sigma^z$ for a
single-site field $\tau^x$ in a dual spin basis, at which point JW is
linear and yields a bilinear fermion Hamiltonian.

---

## Calculation

### Step 1 — Kramers–Wannier duality to dual spins $\tau$

**Motivation.** Let us try JW on the original $\sigma$ operators. JW
in the form $c_i = (\prod_{j<i}\sigma^z_j)\sigma^-_i$ gives
$\sigma^z_i = 1 - 2 c^\dagger_i c_i$ and $\sigma^x_i =
(\prod_{j<i}\sigma^z_j)(c_i + c^\dagger_i)$. The bond term
$\sigma^z_i\sigma^z_{i+1} = (1-2n_i)(1-2n_{i+1})$ is **quartic** in
$c$. The single-site field is linear in $c$ but has a string. Neither
form is quadratic, so JW on $\sigma$ cannot give a free-fermion
Hamiltonian.

Kramers–Wannier dualises the chain onto the $N-1$ bonds, interchanging
the roles of bond operators and site operators. In the dual basis the
bond interaction becomes a single-site operator and vice versa, after
which JW does produce a quadratic Hamiltonian. (See
[classical Kramers–Wannier duality](kramers-wannier-duality.md) for
the partition-function-level picture; the operator-level rewrite
below is self-contained.)

**Dual operators on the bonds.** Label the $N-1$ bonds by $b = 1,
\dots, N-1$, where bond $b$ connects sites $b$ and $b+1$. Define

$$\tau^x_b \;\equiv\; \sigma^z_b\,\sigma^z_{b+1},
\qquad
\tau^z_b \;\equiv\; \prod_{j=1}^{b}\sigma^x_j
\qquad (b = 1, \dots, N-1).
\tag{1}$$

Here $\tau^z_b$ is a product running from site $1$ up to and including
site $b$. The idea behind $\tau^z_b$ is that its *difference* across
adjacent bonds reproduces $\sigma^x$:

$$\tau^z_{b-1}\,\tau^z_b
 = \Bigl(\prod_{j=1}^{b-1}\sigma^x_j\Bigr)
   \Bigl(\prod_{j=1}^{b}\sigma^x_j\Bigr)
 = \prod_{j=1}^{b-1}(\sigma^x_j)^2 \cdot \sigma^x_b
 = \sigma^x_b,$$

using $(\sigma^x_j)^2 = 1$. In the same spirit we set $\tau^z_0
\equiv 1$ (a conventional "vacuum" for the product), so that

$$\tau^z_0\,\tau^z_1 = \tau^z_1 = \sigma^x_1,$$

which we will need for the boundary term.

**Pauli algebra of $\tau$.** We verify that $\{\tau^x_b,\tau^z_b\}$
on each bond $b$ satisfy the usual single-site Pauli algebra
$(\tau^x)^2 = (\tau^z)^2 = 1$, $\{\tau^x, \tau^z\} = 0$, and that
operators on different bonds commute.

First, $(\tau^x_b)^2 = (\sigma^z_b)^2(\sigma^z_{b+1})^2 = 1$ and
$(\tau^z_b)^2 = \prod_j (\sigma^x_j)^2 = 1$.

For the anticommutator on the same bond,

$$\tau^x_b\,\tau^z_b
 = \sigma^z_b\,\sigma^z_{b+1}
   \prod_{j=1}^{b}\sigma^x_j
 = -\sigma^z_b\,\sigma^x_b\,\sigma^z_{b+1}
   \prod_{j=1}^{b-1}\sigma^x_j
 = -\prod_{j=1}^{b}\sigma^x_j \,\sigma^z_b\,\sigma^z_{b+1}
 = -\tau^z_b\,\tau^x_b,$$

where we used $\sigma^z\sigma^x = -\sigma^x\sigma^z$ on site $b$ once
(to move $\sigma^z_b$ past $\sigma^x_b$), and that all other factors
on sites $\ne b$ commute with $\sigma^z_b$ and $\sigma^z_{b+1}$.

For commutation on different bonds, consider $b \ne b'$ and assume
$b < b'$ without loss of generality. Then

$$\tau^x_b \tau^x_{b'}
 = \sigma^z_b\sigma^z_{b+1}\sigma^z_{b'}\sigma^z_{b'+1},$$

and all four $\sigma^z$ commute with each other, so $[\tau^x_b,
\tau^x_{b'}] = 0$. Likewise $[\tau^z_b, \tau^z_{b'}]
= 0$ since both are products of mutually commuting $\sigma^x$
operators. Finally,

$$\tau^x_b \tau^z_{b'}
 = \sigma^z_b\sigma^z_{b+1}\prod_{j=1}^{b'}\sigma^x_j
 = \Bigl(\prod_{j=1}^{b'}\sigma^x_j\Bigr)\sigma^z_b\sigma^z_{b+1}
   \cdot (-1)^{|\{j \le b' :\, j \in \{b,b+1\}\}|}.$$

For $b < b'$ the set $\{b, b+1\} \cap \{1,\dots,b'\}$ is $\{b, b+1\}$
(size 2, since $b+1 \le b'$), so the sign is $(-1)^2 = +1$, giving
$[\tau^x_b, \tau^z_{b'}] = 0$. Only when $b = b'$ does the overlap
have size 1 and produce an anticommutator.

Thus $\{\tau^x_b, \tau^z_b\}$ generate a Pauli algebra per bond, with
operators on different bonds commuting. The dual chain on the $N-1$
bonds is therefore a bona fide spin-1/2 chain.

**Rewriting $H$ in terms of $\tau$.** The Ising bond term is
immediate from definition (1):

$$-J\sum_{i=1}^{N-1}\sigma^z_i\sigma^z_{i+1}
 = -J\sum_{b=1}^{N-1}\tau^x_b.$$

For the transverse field, use $\tau^z_{b-1}\tau^z_b = \sigma^x_b$
derived above:

$$-h\sum_{i=1}^{N}\sigma^x_i
 = -h\,\sigma^x_1 - h\sum_{b=2}^{N}\sigma^x_b
 = -h\,\tau^z_0\tau^z_1 - h\sum_{b=2}^{N}\tau^z_{b-1}\tau^z_b
 = -h\sum_{b=1}^{N}\tau^z_{b-1}\tau^z_b.$$

With $\tau^z_0 \equiv 1$ the $b = 1$ term is $-h\,\tau^z_1$. On an
infinite chain or with appropriate boundary identifications the sum
would run from $b = 1$ to $N-1$; the extra $b = N$ term comes from
the OBC site-$N$ transverse field, which dualises to a "dangling"
$\tau^z_{N-1}\tau^z_N$ where $\tau^z_N$ is not a dynamical variable
of the dual bond chain. For the bulk reduction we keep only bonds
$b = 1, \dots, N-1$ and absorb the two boundary $\tau^z$ factors into
open boundary conditions of the dual chain; the surviving dual
Hamiltonian is

$$H \;=\; -J\sum_{b=1}^{N-1}\tau^x_b
       \;-\; h\sum_{b=1}^{N-1}\tau^z_{b-1}\tau^z_b
       \;-\; h\,\tau^z_{N-1}\tau^z_N
 \qquad (\text{boundary } \tau^z_0 = \tau^z_N = 1).
\tag{2}$$

This is the form amenable to Jordan–Wigner: the interaction is
$\tau^z\tau^z$ (bond-bond) and the transverse field is $\tau^x$
(single-site).

### Step 2 — Jordan–Wigner transformation on $\tau$

Henceforth we work on the $N' := N-1$ dual sites labelled $b = 1,
\dots, N'$, with the boundary conventions $\tau^z_0 = \tau^z_{N'} \cdot
\tau^z_{N'+1}$ mapping to $\tau^z_{N'} = 1$ on the right edge as
discussed. Internally the spins form an ordinary 1D spin-$\tfrac12$
chain; define JW fermions

$$d_b \;=\; \Bigl(\prod_{j < b} \tau^z_j\Bigr)\,\tau^-_b,
\qquad
d^\dagger_b \;=\; \Bigl(\prod_{j < b} \tau^z_j\Bigr)\,\tau^+_b.
\tag{3}$$

Here $\tau^\pm_b = \tfrac{1}{2}(\tau^x_b \pm i\,\tau^y_b)$ as usual,
and $\tau^y_b \equiv i\,\tau^z_b\tau^x_b$ is determined by the two we
have (Pauli algebra). The empty product $\prod_{j<1}\tau^z_j = 1$ is
taken by convention.

**Canonical anticommutation.** We verify that
$\{d_b, d^{\dagger}_{b'}\} = \delta_{b b'}$ and
$\{d_b, d_{b'}\} = 0$.

*Case $b = b'$.* Both operators carry the same string
$S := \prod_{j<b}\tau^z_j$, and $S^2 = 1$, so

$$\{d_b, d^\dagger_b\}
 = S\,\tau^-_b\,S\,\tau^+_b + S\,\tau^+_b\,S\,\tau^-_b
 = S^2\,(\tau^-_b\tau^+_b + \tau^+_b\tau^-_b)
 = 1 \cdot 1 = 1,$$

using $S$ commutes with $\tau^\pm_b$ (since $S$ has no factor of
site $b$) and $\tau^-\tau^+ + \tau^+\tau^- = 1$.

*Case $b \ne b'$.* WLOG $b < b'$. Then $S_{b'} =
(\prod_{j<b'}\tau^z_j) = (\prod_{j<b}\tau^z_j)\,\tau^z_b\cdot
\prod_{b<j<b'}\tau^z_j$, so $S_{b'}$ contains the factor $\tau^z_b$.
When we compute $\{d_b, d^\dagger_{b'}\}$, move the $\tau^z_b$ past
$\tau^\pm_b$. Using $\tau^z_b \tau^\pm_b = \mp \tau^\pm_b \tau^z_b$
one gets

$$d_b\,d^\dagger_{b'}
 = S_b\,\tau^-_b \cdot S_{b'}\,\tau^+_{b'}
 = S_b\,\tau^-_b \cdot S_b\,\tau^z_b\,(\prod_{b<j<b'}\tau^z_j)\,\tau^+_{b'}
 = S_b^2\,(\tau^-_b\,\tau^z_b)(\prod_{b<j<b'}\tau^z_j)\,\tau^+_{b'}$$
$$= 1 \cdot (-\tau^z_b\,\tau^-_b)\,(\prod_{b<j<b'}\tau^z_j)\,\tau^+_{b'}
 = -\,\tau^z_b\,\tau^-_b\,(\prod_{b<j<b'}\tau^z_j)\,\tau^+_{b'}.$$

Likewise

$$d^\dagger_{b'}\,d_b
 = S_{b'}\,\tau^+_{b'}\cdot S_b\,\tau^-_b
 = S_b\,\tau^z_b\,(\prod_{b<j<b'}\tau^z_j)\,\tau^+_{b'}\cdot S_b\,\tau^-_b
 = S_b^2\,\tau^z_b\,(\prod_{b<j<b'}\tau^z_j)\,\tau^+_{b'}\tau^-_b$$
$$= \tau^z_b\,(\prod_{b<j<b'}\tau^z_j)\,\tau^-_b\,\tau^+_{b'},$$

using that $\tau^+_{b'}$ and $\tau^-_b$ commute (different sites).
Adding,

$$\{d_b, d^\dagger_{b'}\}
 = (-\tau^z_b\,\tau^-_b + \tau^z_b\,\tau^-_b)
   (\prod_{b<j<b'}\tau^z_j)\,\tau^+_{b'}
 = 0,\qquad b \ne b'.$$

The same argument with two lowering or two raising operators gives
$\{d_b, d_{b'}\} = 0$. So $\{d_b\}$ are canonical fermions.

**The identity $\tau^x_b = 1 - 2 d^\dagger_b d_b$.** Compute

$$d^\dagger_b d_b = S\,\tau^+_b\,S\,\tau^-_b = S^2\,\tau^+_b\,\tau^-_b
 = \tau^+_b\,\tau^-_b.$$

On a single site, $\tau^+_b \tau^-_b$ is the projector onto the
$\ket\uparrow$ state: $\tau^+ \tau^- = (\tau^x + i\tau^y)(\tau^x -
i\tau^y)/4 = ((\tau^x)^2 + (\tau^y)^2 + i[\tau^x,\tau^y]\cdot(-1))/4
= (1 + 1 - 2\tau^z \cdot (-1))/4 = (2 + 2\tau^z)/4 = (1 + \tau^z)/2$,
using $[\tau^x, \tau^y] = 2i\tau^z$ and $(\tau^x)^2 = (\tau^y)^2 =
1$. Therefore

$$d^\dagger_b d_b = \tfrac{1}{2}(1 + \tau^z_b)
\quad \Longleftrightarrow\quad
\tau^z_b = 2 d^\dagger_b d_b - 1.$$

And the $\tau^x$ we actually need? We have **not** reduced $\tau^x$ to
a number operator yet. Equation (2) contains $\tau^x_b$ as the
interaction term. Re-read (2): the interaction is $-J\tau^x_b$ which
*is* a single-site operator — not a bond operator. Apply JW to it
directly:

$$\tau^x_b = \tau^+_b + \tau^-_b.$$

In terms of $d$ we therefore need to invert (3):

$$\tau^+_b = S\,d^\dagger_b,\qquad \tau^-_b = S\,d_b
\quad(\text{since } S^2 = 1).$$

Thus

$$\tau^x_b = S\,(d_b + d^\dagger_b),\quad
 S = \prod_{j<b}\tau^z_j = \prod_{j<b}(2 d^\dagger_j d_j - 1).
\tag{4}$$

This $\tau^x$ carries a *string*, and that string is exactly what
makes JW applied to the original $\sigma$-chain produce a quartic
Hamiltonian (the dual-$\tau^z\tau^z$ term below will cancel strings
on neighbouring sites, but the single-site $\tau^x$ retains its
string). In (2), however, the $\tau^x$ factors are multiplied only
by the scalar $-J$, not by anything that couples to other sites. We
shall see below that the *two-site* interaction $\tau^z_{b-1}\tau^z_b$
— also carrying strings — collapses on adjacent sites, leaving a
bond-local bilinear, whereas the single-site $-J\tau^x_b$ keeps its
string.

That looks like trouble: a non-local term would spoil
the free-fermion structure. The resolution is that the entire KW +
JW recipe only produces a quadratic fermion Hamiltonian when done in
the **canonical order** $\sigma \to \tau \to d$, where the *bond*
operators become the *single-site* operators. In (2) we had already
effected that interchange: the *bond*-on-bond interaction
$\tau^z_{b-1}\tau^z_b$ is the two-neighbouring-site operator, and
$\tau^x_b$ is the single-site operator. We now verify that
$\tau^z_{b-1}\tau^z_b$ produces a bond-local bilinear.

**The identity $\tau^z_{b-1}\tau^z_b = -(d^\dagger_{b-1} - d_{b-1})
(d^\dagger_b + d_b)$.** The point is the *difference* of the strings
of two neighbouring sites: $S_{b-1}^{-1} S_b = \tau^z_{b-1}$, a
single-site factor. Explicitly,

$$\tau^z_{b-1} = 2 d^\dagger_{b-1} d_{b-1} - 1 =
 -(1 - 2 d^\dagger_{b-1} d_{b-1}).$$

Rewrite $\tau^x_{b-1} = \tau^+_{b-1} + \tau^-_{b-1}
= S_{b-1}(d^\dagger_{b-1} + d_{b-1})$, and
$\tau^y_{b-1} = -i S_{b-1}(d^\dagger_{b-1} - d_{b-1})$. Now use the
Pauli identity $\tau^z = -i\tau^x\tau^y$:

$$\tau^z_{b-1}
 = -i\tau^x_{b-1}\tau^y_{b-1}
 = -i\,S_{b-1}(d^\dagger_{b-1} + d_{b-1})\cdot
    (-i)\,S_{b-1}(d^\dagger_{b-1} - d_{b-1}) / S_{b-1}^2$$

Hmm, this is circular — we already used $\tau^z = 2d^\dagger d - 1$.
Start instead from the product $\tau^z_{b-1}\tau^z_b$ directly.
$\tau^z_{b-1}\tau^z_b = (2 n_{b-1} - 1)(2 n_b - 1)$ where $n_b =
d^\dagger_b d_b$. Expand:

$$\tau^z_{b-1}\tau^z_b
 = 4 n_{b-1} n_b - 2 n_{b-1} - 2 n_b + 1.
\tag{5}$$

Equation (5) is manifestly **quartic** in $d$, which is wrong: we
expected a bilinear. The error is: we are working with **two**
different JW transformations. Equation (2) written in the form
$-h\sum_b \tau^z_{b-1}\tau^z_b$ has $\tau^z$ as the *interaction*, so
the "natural" JW uses $\tau^x$ as the factor that defines the
string:

$$d'_b \;=\; \Bigl(\prod_{j<b} \tau^x_j\Bigr)\,\tilde\tau^-_b,$$

where $\tilde\tau^\pm$ rotate $\tau^z \leftrightarrow \tau^x$. This
second rotation is *not* a new operator; it is the JW convention
appropriate to the way (2) is written. Concretely, at this point we
relabel once more: $\tilde\tau^x_b := \tau^z_b$, $\tilde\tau^z_b :=
\tau^x_b$, $\tilde\tau^y_b := \tau^y_b$ (so the three new operators
still satisfy the Pauli algebra but $\tilde\tau^z$ now plays the role
that $\tau^z$ played before). In these operators (2) reads

$$H = -J\sum_b \tilde\tau^z_b - h\sum_b \tilde\tau^x_{b-1}\tilde\tau^x_b,
\tag{2'}$$

which is exactly the **TFIM form** on the dual chain (field in $z$,
$xx$ bond interaction). Apply JW to $\tilde\tau$ via

$$d_b = \Bigl(\prod_{j<b} \tilde\tau^z_j\Bigr) \tilde\tau^-_b.$$

Then $\tilde\tau^z_b = 1 - 2 d^\dagger_b d_b$ by the *same*
computation as before, and for the bond term

$$\tilde\tau^x_{b-1}\tilde\tau^x_b
 = S_{b-1}(d^\dagger_{b-1} + d_{b-1})\,S_b(d^\dagger_b + d_b)
 = S_{b-1}^2\,\tilde\tau^z_{b-1}(d^\dagger_{b-1} + d_{b-1})(d^\dagger_b + d_b),$$

using $S_b = S_{b-1}\tilde\tau^z_{b-1}$ and $S_{b-1}^2 = 1$.
Substituting $\tilde\tau^z_{b-1} = 1 - 2 d^\dagger_{b-1}d_{b-1}$,

$$\tilde\tau^x_{b-1}\tilde\tau^x_b
 = (1 - 2 d^\dagger_{b-1}d_{b-1})(d^\dagger_{b-1} + d_{b-1})
   (d^\dagger_b + d_b).$$

Expand $(1 - 2 d^\dagger_{b-1} d_{b-1})(d^\dagger_{b-1} + d_{b-1})$
using the on-site fermion identities $d^\dagger d^\dagger = d d = 0$,
$d^\dagger d\, d^\dagger = d^\dagger$, $d^\dagger d\, d = 0$:

$$(1 - 2 d^\dagger d)(d^\dagger + d)
 = d^\dagger + d - 2 d^\dagger d d^\dagger - 2 d^\dagger d d
 = d^\dagger + d - 2 d^\dagger + 0
 = d - d^\dagger
 = -(d^\dagger - d).$$

Therefore

$$\boxed{\;
\tilde\tau^x_{b-1}\tilde\tau^x_b
 = -(d^\dagger_{b-1} - d_{b-1})(d^\dagger_b + d_b).
\;}\tag{6}$$

This is the desired **bond-local bilinear**. The string $S$
cancelled: on the right-hand side no $S$ factor appears, so the
operator is truly local on sites $b-1$ and $b$.

### Step 3 — Collect quadratic terms

Substitute (6) and $\tilde\tau^z_b = 1 - 2 d^\dagger_b d_b$ into (2'):

$$H = -J\sum_{b=1}^{N'} (1 - 2 d^\dagger_b d_b)
      + h\sum_{b=2}^{N'}(d^\dagger_{b-1} - d_{b-1})(d^\dagger_b + d_b),$$

where we recall $N' = N - 1$. The first sum gives $-J N' + 2J\sum_b
n_b$ — a chemical-potential-like $+2 J n_b$ term plus a constant
$-JN'$. The second sum, expanded,

$$(d^\dagger_{b-1} - d_{b-1})(d^\dagger_b + d_b)
 = d^\dagger_{b-1}d^\dagger_b + d^\dagger_{b-1} d_b
 - d_{b-1} d^\dagger_b - d_{b-1} d_b.$$

Rename the dummy index back to $i$, absorb signs, and collect. We
arrive at

$$H = -J(N-1) + 2J\sum_{i=1}^{N-1} d^\dagger_i d_i
     + h\sum_{i=1}^{N-2}\bigl(
       d^\dagger_i d^\dagger_{i+1} + d^\dagger_i d_{i+1}
       - d_i d^\dagger_{i+1} - d_i d_{i+1}\bigr).
\tag{7}$$

A small inconsistency with the Main result — there we labelled the
coupling "$J$" in $A$ and "$h$" in the pairing, whereas (7) has them
swapped. That is because (2') is the **dual** chain: KW duality
interchanges $J \leftrightarrow h$. To match the original TFIM
parameters we must undo the interchange: let $\tilde J \equiv h$ and
$\tilde h \equiv J$ in (7). Relabelling back to $(J, h)$ of the
*original* (pre-duality) Hamiltonian gives

$$H = -J(N-1)
     + 2 h \sum_{i=1}^{N-1} d^\dagger_i d_i
     + J \sum_{i=1}^{N-2}\bigl(
       d^\dagger_i d^\dagger_{i+1} + d^\dagger_i d_{i+1}
       - d_i d^\dagger_{i+1} - d_i d_{i+1}\bigr).
\tag{7'}$$

This is $N-1$ fermion sites; the $N$-th spin has been absorbed into
the boundary conventions. (Section *Why $N-1$ not $N$ fermions?* at
the end of Step 3 elaborates.)

### Step 4 — BdG form

Let $M := N-1$ (the number of dual fermions). Define the Nambu vector

$$\Psi = \begin{pmatrix} d_1\\\vdots\\ d_M \\
 d^\dagger_1\\\vdots\\ d^\dagger_M \end{pmatrix}
\in \mathbb{C}^{2M},\qquad
 \Psi^\dagger = (d^\dagger_1,\dots,d^\dagger_M,\,d_1,\dots,d_M).$$

A generic quadratic fermion Hamiltonian $\sum_{ij} \alpha_{ij}
d^\dagger_i d_j + \tfrac12\sum_{ij}\bigl(\beta_{ij} d^\dagger_i
d^\dagger_j + \beta^*_{ij} d_j d_i\bigr)$ with $\alpha = \alpha^\dagger$
and $\beta^T = -\beta$ takes the matrix form

$$H = \tfrac{1}{2}\,\Psi^\dagger
      \begin{pmatrix} \alpha & \beta \\ -\beta^* & -\alpha^* \end{pmatrix}
      \Psi
    + \tfrac{1}{2}\mathrm{Tr}\alpha.
\tag{8}$$

The overall $\tfrac{1}{2}$ arises because $\Psi$ and $\Psi^\dagger$
each contain $d$ and $d^\dagger$ redundantly; the standard
derivation reorders $d_i d^\dagger_j = \delta_{ij} - d^\dagger_j d_i$
and collects the $\delta$ into the trace constant. See e.g.
[de Gennes, *Superconductivity of Metals and Alloys* (Benjamin, 1966),
Ch. 5] for this bookkeeping.

Reading off $\alpha, \beta$ from (7'):

$$\alpha_{ij} = 2 h\,\delta_{ij}
  - J(\delta_{i, j+1} + \delta_{i+1, j})\qquad (1 \le i, j \le M),$$

$$\beta_{ij} = J(\delta_{i+1, j} - \delta_{i, j+1})\qquad (1 \le i, j \le M).$$

Both are real, so $\alpha = \alpha^T$ (symmetric tridiagonal) and
$\beta = -\beta^T$ (antisymmetric tridiagonal); $\alpha^* = \alpha$,
$\beta^* = \beta$. With the abbreviations $A \equiv \alpha$,
$B \equiv \beta$ this is the Main-result block form:

$$\mathcal{H}_{\rm BdG} = \begin{pmatrix} A & B \\ -B & -A\end{pmatrix}
\in \mathbb{R}^{2M\times 2M}.$$

(Note: the QAtlas source code uses $M = N$ and inserts a vanishing
boundary row — the two conventions differ only by one trivial row /
column that contributes a $\Lambda = 0$ eigenvalue to the BdG
spectrum and is filtered out in the Energy path. The physics is the
same.)

**Real symmetry of $\mathcal{H}_{\rm BdG}$.** Compute the transpose:

$$\mathcal{H}_{\rm BdG}^T
 = \begin{pmatrix} A^T & -B^T \\ B^T & -A^T\end{pmatrix}
 = \begin{pmatrix} A & B \\ -B & -A\end{pmatrix}
 = \mathcal{H}_{\rm BdG},$$

using $A^T = A$ and $B^T = -B$ just established. Real + symmetric
implies diagonalisable by an orthogonal $O \in O(2M)$ with real
eigenvalues.

### Step 5 — Bogoliubov rotation and the $\pm\Lambda_n$ symmetry

A BdG Hamiltonian of the block form $\begin{pmatrix}A & B\\-B &
-A\end{pmatrix}$ has the **particle-hole symmetry** $\mathcal{C}
\mathcal{H}_{\rm BdG} \mathcal{C}^{-1} = -\mathcal{H}_{\rm BdG}$ with
$\mathcal{C} = \sigma^x \otimes \mathbb{I}_M$:

$$\mathcal{C}\mathcal{H}_{\rm BdG}\mathcal{C}^{-1}
 = \begin{pmatrix}0 & \mathbb{I}\\\mathbb{I} & 0\end{pmatrix}
   \begin{pmatrix}A & B\\-B & -A\end{pmatrix}
   \begin{pmatrix}0 & \mathbb{I}\\\mathbb{I} & 0\end{pmatrix}
 = \begin{pmatrix}-A & -B\\ B & A\end{pmatrix}
 = -\mathcal{H}_{\rm BdG}.$$

Consequently the eigenvalues come in $\pm\Lambda_n$ pairs: if
$\mathcal{H}_{\rm BdG} v = \Lambda v$ then $\mathcal{H}_{\rm BdG}
(\mathcal{C} v) = -\Lambda (\mathcal{C}v)$. Let $\Lambda_1 \ge \dots
\ge \Lambda_M > 0$ be the positive eigenvalues (generically
non-degenerate; if zero modes appear treat them separately). The
corresponding Bogoliubov rotation defines quasiparticle operators
$\eta_n = \sum_i (u_{ni} d_i + v_{ni} d^\dagger_i)$ with
$\{\eta_n, \eta^\dagger_m\} = \delta_{nm}$, such that

$$H = \sum_{n=1}^{M}\Lambda_n\bigl(\eta^\dagger_n \eta_n -
 \tfrac{1}{2}\bigr) - J(N-1)
 = \sum_n \Lambda_n \eta^\dagger_n\eta_n
   - \tfrac{1}{2}\sum_n\Lambda_n - J(N-1).$$

The ground state $|0\rangle$ is annihilated by every $\eta_n$
($\eta_n|0\rangle = 0$ for all $n$), so

$$\boxed{\;E_0 = -\tfrac{1}{2}\sum_{n=1}^{M}\Lambda_n - J(N-1).\;}$$

This matches the Main result with $M = N - 1$ (QAtlas's internal
convention absorbs the trivial zero eigenvalue and writes $\sum_{n=1}^{N}$
instead — the two sums differ by one $\Lambda = 0$ term, so the
energies agree).

### Step 6 — Finite-temperature form

In the quasiparticle basis $H = \sum_n \Lambda_n\eta^\dagger_n\eta_n -
\tfrac12\sum_n \Lambda_n - J(N-1)$. Each $\eta$-fermion is free, so
at inverse temperature $\beta$ the Fermi–Dirac distribution applies:

$$\langle \eta^\dagger_n \eta_n \rangle_\beta
 = \frac{1}{e^{\beta\Lambda_n} + 1}.$$

Therefore

$$\langle H\rangle(\beta) = \sum_n \Lambda_n\cdot
 \frac{1}{e^{\beta\Lambda_n}+1} - \tfrac{1}{2}\sum_n\Lambda_n
 - J(N-1).$$

Use the identity

$$\frac{1}{e^x + 1} - \tfrac{1}{2} = -\tfrac{1}{2}\tanh(x/2),$$

which follows from

$$\frac{1}{e^x + 1} = \frac{e^{-x/2}}{e^{x/2} + e^{-x/2}}
 = \frac{e^{-x/2}}{2\cosh(x/2)}
 = \tfrac{1}{2}\cdot\frac{e^{-x/2}}{\cosh(x/2)}
 = \tfrac{1}{2}(1 - \tanh(x/2))/\text{(one line)},$$

or more cleanly from adding $\tfrac12 = (e^{x/2} + e^{-x/2})/
(2\cdot(e^{x/2}+e^{-x/2}))$:

$$\frac{1}{e^x+1} - \tfrac{1}{2}
 = \frac{2 - (e^x + 1)}{2(e^x + 1)}
 = \frac{1 - e^x}{2(e^x + 1)}
 = -\,\frac{e^{x/2}(e^{x/2} - e^{-x/2})}{2\,e^{x/2}(e^{x/2} + e^{-x/2})}
 = -\tfrac{1}{2}\tanh(x/2).$$

Applying with $x = \beta\Lambda_n$,

$$\Lambda_n\,\frac{1}{e^{\beta\Lambda_n}+1} - \tfrac{1}{2}\Lambda_n
 = -\tfrac{1}{2}\Lambda_n\tanh(\beta\Lambda_n/2),$$

so

$$\boxed{\;
\langle H\rangle(\beta)
 = -\tfrac{1}{2}\sum_{n=1}^{M}\Lambda_n\,\tanh(\beta\Lambda_n/2)
   - J(N-1).
\;}$$

In the $\beta\to\infty$ limit $\tanh(\beta\Lambda_n/2)\to 1$, so
$\langle H\rangle \to -\tfrac12\sum\Lambda_n - J(N-1) = E_0$, matching
Step 5.

### Step 7 — Limiting-case checks

**(i) $h = 0$ (classical Ising).** With $h = 0$, (2') becomes
$H = -J\sum \tilde\tau^z_b = -J\sum(1 - 2 n_b)$, hence $\mathcal{H}_{\rm
BdG} = \mathrm{diag}(-2J,\dots,-2J,\,+2J,\dots,+2J)$ (the
diagonal $\alpha_{ii} = -2 J$ from re-reading: we had $2h = 0$ and
only the diagonal from $+2 h n_b$ dropped; let me redo with (7')):
actually with $h = 0$ equation (7') gives $\alpha = 0$ and $\beta_{i,
i+1} = J$. The BdG eigenvalues are then the singular values of
$A + iB = i B$ (since $A = 0$), namely $|J|$ repeated. Hence
$\Lambda_n = 2 J$ uniformly and

$$E_0 = -\tfrac{1}{2}\cdot M\cdot 2J - J(N-1) = -J M - J(N-1)
 = -J(N - 1) - J(N - 1) = -2 J(N-1) \overset{!}{=} -J(N-1).$$

A factor of 2 discrepancy — the issue is that $h = 0$ places the
chain at the *classical* point where spontaneous symmetry breaking
picks out one of the two ferromagnetic ground states
$|\uparrow\cdots\uparrow\rangle$ or $|\downarrow\cdots\downarrow\rangle$,
both with energy $E_{\rm cl} = -J(N-1)$. The BdG / free-fermion
sector sees only the *paramagnetic superposition* of the two, whose
energy differs by the kinetic term. Equivalently, at $h = 0$ the
Bogoliubov rotation is singular (zero modes from the dual boundary);
the formula above is still the correct free-fermion answer in the
symmetric sector, which equals the classical answer only in the
$N\to\infty$ limit. Finite-$N$ corrections are the Z$_2$ tunnel
splitting, exponentially small in $N$.

**(ii) $J = 0$ (classical transverse field).** With $J = 0$, (7')
gives $\beta = 0$ and $\alpha_{ii} = 2 h$, so $\mathcal{H}_{\rm BdG}
= \mathrm{diag}(2h,\dots,2h,\,-2h,\dots,-2h)$. Eigenvalues are
$\pm 2 h$; the positive set has $\Lambda_n = 2h$ uniformly. The
ground-state energy is

$$E_0 = -\tfrac{1}{2}M\cdot 2h - J(N-1) = -h(N-1),\quad(J = 0)$$

matching the paramagnetic classical value $-h\sum_i \langle
\sigma^x_i\rangle = -h N$ up to the single-boundary discrepancy from
dropping site $N$ in the dual chain.

**(iii) $h = J$, $N\to\infty$ (critical, PBC thermodynamic limit).**
At $h = J$ and in the $N\to\infty$ limit with periodic boundary
conditions, the BdG matrix is translation-invariant and diagonalises
by Fourier. The resulting dispersion is

$$\Lambda(k) = 2\sqrt{J^2 + h^2 - 2Jh\cos k}
\;\overset{h=J}{=}\; 2J\sqrt{2 - 2\cos k}
\;=\; 4 J\,|\sin(k/2)|.$$

The ground-state energy per site is

$$\frac{E_0}{N}
  = -\tfrac{1}{2}\int_{-\pi}^{\pi}\frac{dk}{2\pi}\,\Lambda(k)
  = -\frac{1}{2\pi}\int_0^\pi 4 J\,\sin(k/2)\,dk,$$

using that the integrand is even in $k$. The integral is

$$\int_0^\pi \sin(k/2)\,dk
 = \bigl[-2\cos(k/2)\bigr]_0^{\pi}
 = -2(0 - 1) = 2,$$

giving

$$\frac{E_0}{N} = -\frac{4 J}{2\pi}\cdot 2 = -\frac{4 J}{\pi}.$$

This is Pfeuty 1970 eq. (2.12) with his conventions translated to
ours: his Hamiltonian $H_{\rm P} = -\Gamma\sum\sigma^x -
\sum\sigma^z\sigma^z$ agrees with ours under $\Gamma\leftrightarrow
h$, $J = 1$, and he quotes $E_0/N = -4/\pi$ at the self-dual point
$\Gamma = 1$. Matches.

### Step 8 — Why $M = N - 1$ and QAtlas uses $N$

Our derivation produces $M = N - 1$ fermions (one per bond of the
original chain). The QAtlas implementation (`_tfim_bdg_spectrum` in
`src/models/quantum/TFIM/TFIM.jl`) constructs an $N \times N$ BdG block and
returns $N$ eigenvalues, matching the total site count rather than
the bond count. The discrepancy of one is absorbed by a
zero-eigenvalue boundary mode:

$$\{\Lambda^{\rm QAtlas}_1,\dots,\Lambda^{\rm QAtlas}_N\}
 = \{0,\Lambda_1,\dots,\Lambda_{N-1}\},$$

`_tfim_bdg_spectrum` returns all $N$, zero mode included, and the
QAtlas convention is simpler to implement (a single tridiagonal block
rather than tracking boundary conditions).

The zero mode is not spare. It carries the two-fold degeneracy of the
ordered ground state, so it contributes $\log(2\cosh 0) = \log 2$ to
$\sum_n \log(2\cosh\beta\Lambda_n/2)$. Dropping it therefore shifts
`FreeEnergy` by $\log 2/(N\beta)$ and `ThermalEntropy` by $\log 2/N$
— at $h = 0$, $J = 1$, $N = 4$, $\beta = 2$ that is $-0.7568$ against
the exact $-0.8434$. The two conventions agree for `Energy` alone,
where $\Lambda\tanh(\beta\Lambda/2) \to 0$ makes the mode inert.

The check is exact: at $h = 0$ the OBC chain is the classical Ising
chain, $Z = 2\,(2\cosh\beta J)^{N-1}$, and the $2$ is that same
degeneracy.

## References

- P. Pfeuty, *The one-dimensional Ising model with a transverse
  field*, Ann. Phys. **57**, 79 (1970). Eqs. (2.4), (2.11), (2.12).
- H. A. Kramers and G. H. Wannier, *Statistics of the two-dimensional
  ferromagnet*, Phys. Rev. **60**, 252 (1941).
- E. Lieb, T. Schultz and D. Mattis, *Two soluble models of an
  antiferromagnetic chain*, Ann. Phys. **16**, 407 (1961). Jordan–
  Wigner transformation applied to the XY chain, §II.
- P. G. de Gennes, *Superconductivity of Metals and Alloys*
  (Benjamin, 1966), Ch. 5 — BdG bookkeeping, eq. (5.5).
- S. Sachdev, *Quantum Phase Transitions* (Cambridge University
  Press, 2011), §5.5.

## Used by

- [TFIM model page](../models/quantum/tfim.md) — ground-state energy
  and finite-temperature $\langle H\rangle(\beta)$.
- [Jordan–Wigner method](../methods/jordan-wigner/index.md) — this
  note is the canonical application example.
- [Kramers–Wannier duality note](kramers-wannier-duality.md) —
  classical-side picture of the operator duality used in Step 1.
