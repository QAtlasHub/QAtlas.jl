# CFT: Virasoro Minimal Models and WZW SU(2)$_k$

QAtlas exposes two parametric families of two-dimensional rational
conformal field theories as concrete dispatch tags:

- [`MinimalModel(p, p_prime)`](@ref) — unitary and non-unitary
  Virasoro minimal models $\mathcal{M}(p, p^\prime)$.
- [`WZWSU2(k)`](@ref) — Wess–Zumino–Witten models with affine
  $\widehat{\mathfrak{su}}(2)_k$ symmetry.

Both expose their CFT data as exact `Rational{Int}` values via
`fetch` on [`CentralCharge`](@ref),
[`ConformalWeights`](@ref), and (minimal models only)
[`PrimaryFields`](@ref).

---

## Virasoro Minimal Models $\mathcal{M}(p, p^\prime)$

Following Belavin, Polyakov, Zamolodchikov (BPZ, 1984) the central
charge of the minimal model labelled by coprime integers
$p > p^\prime \geq 2$ is

$$c(p, p^\prime) = 1 - \frac{6 (p - p^\prime)^2}{p \, p^\prime}.$$

Primary fields are labelled by $(r, s)$ with $1 \leq r \leq p^\prime - 1$,
$1 \leq s \leq p - 1$, with conformal weight given by the **Kac formula**:

$$h_{r, s}(p, p^\prime) = \frac{(p \, r - p^\prime s)^2 - (p - p^\prime)^2}{4 \, p \, p^\prime}.$$

The Kac table has the symmetry $h_{r, s} = h_{p^\prime - r, p - s}$,
so the number of distinct primaries is $(p - 1)(p^\prime - 1) / 2$.

### Special cases

| Model                    | $(p, p^\prime)$ | $c$       | Primary weights $h$            |
|:-------------------------|:----------------|:----------|:--------------------------------|
| Yang–Lee (non-unitary)   | $(5, 2)$        | $-22/5$   | $0, -1/5$                       |
| **Ising** $\mathcal{M}(4,3)$ | $(4, 3)$    | $1/2$     | $h_{1,1}=0,\ h_{1,2}=1/16,\ h_{2,1}=1/2$ |
| Tricritical Ising        | $(5, 4)$        | $7/10$    | $0,\ 3/80,\ 1/10,\ 7/16,\ 3/5,\ 3/2$ |
| 3-state Potts (chiral)   | $(6, 5)$        | $4/5$     | 10 primaries                    |

### Usage

```julia
using QAtlas

m = MinimalModel(4, 3)
fetch(m, CentralCharge())                  # 1//2
fetch(m, ConformalWeights(); r=1, s=2)     # 1//16  (σ)
fetch(m, ConformalWeights(); r=2, s=1)     # 1//2   (ε)
fetch(m, PrimaryFields())                  # 3 NamedTuples (r, s, h)
```

### Cross-check with `Universality(:Ising)`

The `c = 1//2` stored on `Universality(:Ising)`'s exponent
NamedTuple at $d = 2$ agrees identically with
`fetch(MinimalModel(4, 3), CentralCharge())`:

```julia
fetch(MinimalModel(4, 3), CentralCharge()) ==
    QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).c
# true
```

This is a sanity check that the CFT-side parametric data and the
universality-class table agree as `Rational{Int}` values.

### Validation

`MinimalModel(p, p_prime)` validates its arguments at construction
time and throws `DomainError` for any of:

- $p^\prime < 2$,
- $p \leq p^\prime$,
- $\gcd(p, p^\prime) \neq 1$ (non-coprime).

`fetch(::MinimalModel, ConformalWeights(); r, s)` throws
`DomainError` for $(r, s)$ outside the fundamental rectangle
$1 \leq r \leq p^\prime - 1$, $1 \leq s \leq p - 1$.

---

## WZW SU(2)$_k$

Following Knizhnik–Zamolodchikov (1984) the Sugawara central charge
of the level-$k$ WZW model on the affine algebra
$\widehat{\mathfrak{su}}(2)_k$ is

$$c(k) = \frac{3 k}{k + 2}, \qquad k = 1, 2, 3, \dots$$

The primary fields are labelled by SU(2) spin
$j \in \{0, 1/2, 1, 3/2, \dots, k/2\}$ with conformal weight

$$h_j = \frac{j (j + 1)}{k + 2}.$$

### Special cases

| Level | $c$  | Spectrum $j$         | Notes                                                |
|:-----:|:----:|:---------------------|:-----------------------------------------------------|
| 1     | 1    | $0, 1/2$             | Free boson at SU(2)-symmetric radius; low-energy theory of the spin-$1/2$ Heisenberg AFM (Affleck 1989) |
| 2     | 3/2  | $0, 1/2, 1$          | 3 free Majorana fermions (smallest $\mathcal{N}{=}1$ super-Virasoro minimal model) |
| 3     | 9/5  | $0, 1/2, 1, 3/2$     |                                                      |

### Usage

```julia
fetch(WZWSU2(1), CentralCharge())                    # 1//1
fetch(WZWSU2(1), ConformalWeights(); j=1//2)         # 1//4
fetch(WZWSU2(2), ConformalWeights(); j=1)            # 1//2
```

### Validation

`WZWSU2(k)` requires $k \geq 1$ and throws `DomainError` otherwise.

`fetch(::WZWSU2, ConformalWeights(); j)` requires `j` to be a
non-negative half-integer (`Integer` or `Rational{Int}` with
$2 j \in \mathbb{Z}_{\geq 0}$) with $0 \leq j \leq k/2$.  Floats and
non-half-integer rationals (e.g. `1//3`) raise `DomainError`.

---

## References

- A. A. Belavin, A. M. Polyakov, A. B. Zamolodchikov,
  *Infinite conformal symmetry in two-dimensional quantum field theory*,
  Nucl. Phys. **B 241**, 333 (1984).
- D. Friedan, Z. Qiu, S. Shenker,
  *Conformal invariance, unitarity, and critical exponents in two dimensions*,
  Phys. Rev. Lett. **52**, 1575 (1984).
- V. G. Knizhnik, A. B. Zamolodchikov,
  *Current algebra and Wess–Zumino model in two dimensions*,
  Nucl. Phys. **B 247**, 83 (1984).
- E. Witten,
  *Non-abelian bosonization in two dimensions*,
  Comm. Math. Phys. **92**, 455 (1984).
- P. Di Francesco, P. Mathieu, D. Sénéchal,
  *Conformal Field Theory* (Springer, 1997), Ch. 7 (minimal models),
  Ch. 15 (WZW).
- I. Affleck,
  *Quantum spin chains and the Haldane gap*,
  J. Phys.: Condens. Matter **1**, 3047 (1989).
