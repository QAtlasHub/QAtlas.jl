# XX (Δ = 0) Free-Fermion Quench Observables

This page documents the closed-form quench observables exposed by
`QAtlas.fetch(::XXZ1D, ::LoschmidtRateFunction, ::Infinite; ...)`
introduced in issue #148 (phase 1).

## Setup

The XXZ chain in the spin convention `Sᵅ = σᵅ/2`,

```math
H_{\text{XX}}(J) = J \sum_i \bigl( S^x_i S^x_{i+1} + S^y_i S^y_{i+1} \bigr) ,
```

at `Δ = 0` Jordan–Wigner-transforms (after the standard
nearest-neighbour string cancellation) into a tight-binding fermion
chain **without pairing** and at zero chemical potential,

```math
H_{\text{JW}}(J) = \tfrac{J}{2} \sum_i \bigl( c^\dagger_i c_{i+1} + \text{h.c.} \bigr) ,
\qquad
\varepsilon_J(k) = J \cos k .
```

The half-filled ground state is the Slater determinant of all modes
with `ε_J(k) < 0`,

```math
|\mathrm{GS}(J)\rangle
  = \prod_{k\,:\,J \cos k < 0} c^\dagger(k) \, |\varnothing\rangle .
```

Because there is no Bogoliubov pairing, `H_XX(J₀)` and `H_XX(J_f)` are
diagonalised in the **same** plane-wave basis `c(k)` for every choice
of `J₀, J_f`.

## Loschmidt rate function

The Loschmidt rate of the quench `H_XX(J₀) → H_XX(J_f)` from the
initial ground state `|ψ₀⟩ = |GS(J₀)⟩` is, by definition,

```math
\lambda(t)
  = -\lim_{N \to \infty} \frac{1}{N} \log
    \bigl| \langle \psi_0 | e^{-i H_f t} | \psi_0 \rangle \bigr|^2 .
```

In the diagonal basis the amplitude factorises mode by mode:

```math
\langle \psi_0 | e^{-i H_f t} | \psi_0 \rangle
  = \prod_k \langle n_0(k) |
      e^{-i \varepsilon_{J_f}(k) t (\hat n_k - \tfrac12)}
    | n_0(k) \rangle
  = \prod_k \exp\!\left\{
      -i \varepsilon_{J_f}(k) t \cdot \bigl( n_0(k) - \tfrac12 \bigr)
    \right\}
```

with `n_0(k) = Θ(-ε_{J₀}(k))` the initial-state occupation.  The
**modulus is identically 1**, so

```math
\lambda(t) \equiv 0
\qquad
\text{whenever } \operatorname{sgn} J_0 = \operatorname{sgn} J_f .
```

## Fermi-sea topology

The Fermi sea `{k : J cos k < 0}` depends on `sgn J` only.  Three
regimes appear:

| `(sgn J₀, sgn J_f)` | Fermi sea | `λ(t)`               | API behaviour       |
| ------------------- | --------- | -------------------- | ------------------- |
| `(+,+)` or `(-,-)`  | identical | `0` for every `t`    | returns `0.0`        |
| `(+,-)` or `(-,+)`  | complementary | `+∞` (Anderson orth.) | returns `Inf` + `@warn` |
| `(0, 0)`            | flat band, no dynamics | `0`     | returns `0.0`        |
| `(0, ±)` or `(±, 0)` | one flat side, GS degenerate | undefined | returns `NaN` + `@warn` |

The sign-flip case is the well-known orthogonality catastrophe: any
two Slater determinants with complementary occupied-mode sets are
exactly orthogonal in the thermodynamic limit, hence
`|⟨ψ₀ | ψ_f⟩| = 0` and `λ(t) = +∞` for every `t`.

## Why this is degenerate (and what's deferred)

The Calabrese–Essler–Fagotti analysis of XX quench dynamics
(*J. Stat. Mech.* (2012) P07016) treats initial states like the Néel
state or a dimerised state, which are **not Gaussian in the same
fermion basis** as the post-quench Hamiltonian.  Such states induce a
non-trivial single-particle Bogoliubov rotation at the quench
instant, and the Loschmidt amplitude becomes the textbook integral

```math
\lambda(t)
  = -\int_0^\pi \frac{dk}{\pi}
      \log \!\left|
        \cos^2(\Delta\varphi_k) +
        \sin^2(\Delta\varphi_k)\, e^{2 i \varepsilon_{J_f}(k) t}
      \right|^2 .
```

The current `XXZ1D` model carries only `(J, Δ)` with no magnetic
field, dimerisation, or staggered-state machinery, so the only
XX → XX quench expressible at present is `|GS(J₀)⟩ → e^{-iH_f t} |GS(J₀)⟩`,
which is the degenerate `(λ ≡ 0 / +∞)` case derived above.  Phase 2
will lift this restriction by adding either

* a magnetic-field generalisation of `XXZ1D` (so `H₀` and `H_f` differ
  in chemical potential), or
* a separate `XYModel` carrying the pairing γ (so the quench rotates
  the Bogoliubov modes),

at which point the closed-form Calabrese–Essler–Fagotti integral
becomes the meaningful return value.

## API

```julia
fetch(model_f::XXZ1D,
      ::LoschmidtRateFunction,
      ::Infinite;
      initial::XXZ1D,
      t::Real) -> Float64
```

A `Δ ≠ 0` model on either side raises `DomainError`.  The kwarg
`initial` is the initial-state Hamiltonian whose ground state is taken
as `|ψ₀⟩`; `t` is the real evolution time.

## Example

```julia-repl
julia> using QAtlas

julia> m_f = XXZ1D(; J=1.0, Δ=0.0);

julia> m_0 = XXZ1D(; J=0.5, Δ=0.0);

julia> fetch(m_f, LoschmidtRateFunction(), Infinite();
             initial=m_0, t=1.0)
0.0
```

## References

* P. Calabrese, F.H.L. Essler, M. Fagotti,
  *J. Stat. Mech.* (2012) P07016.
* M. Heyl, A. Polkovnikov, S. Kehrein,
  *Phys. Rev. Lett.* 110, 135704 (2013).
* F.H.L. Essler, M. Fagotti,
  *J. Stat. Mech.* (2016) 064002.
