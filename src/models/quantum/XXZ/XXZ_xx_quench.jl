# ─────────────────────────────────────────────────────────────────────────────
# XXZ chain at Δ = 0 (XX / free fermion) — quench observables
# in the thermodynamic limit (Infinite()).
#
# This file is the issue-#148 "phase 1" companion to
# `XXZ_xx_infinite.jl` (PR #157) — both branches add Δ = 0 fetch
# methods that are guarded by `_xx_is_free_fermion`.
#
# ── Setup ────────────────────────────────────────────────────────────────────
#
# In the spin convention `Sᵅ = σᵅ/2` used throughout XXZ.jl, the XX
# chain
#
#     H_XX(J) = J Σᵢ (Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁)
#
# Jordan–Wigner-transforms (after the standard nearest-neighbour string
# cancellation) into a tight-binding fermion chain with **no pairing
# term** and zero chemical potential,
#
#     H_JW(J) = (J/2) Σᵢ (cᵢ† cᵢ₊₁ + h.c.).
#
# Its single-particle dispersion is
#
#     ε_J(k) = J cos(k),     k ∈ [-π, π]                          (★)
#
# at half filling (the µ = 0 line that the spin-1/2 zero-magnetization
# ground state sits on).  The ground state is the Slater determinant
# obtained by filling all modes with `ε_J(k) < 0`,
#
#     |GS(J)⟩ = ∏_{k : J cos k < 0} c†(k) |∅⟩.
#
# Because there is **no Bogoliubov pairing**, single-particle modes do
# not mix between two XX Hamiltonians H_XX(J₀) and H_XX(J_f): the
# operators `c(k)` and `c†(k)` diagonalise both Hamiltonians
# simultaneously, and the only thing that the quench can change is
# **which** modes are below the Fermi level.
#
# ── Loschmidt amplitude under an XX → XX quench ──────────────────────────────
#
# Prepare the system in `|ψ₀⟩ = |GS(J₀)⟩` and evolve under
# `H_f = H_XX(J_f)`.  In the diagonal Bogoliubov-free basis the
# amplitude factorises mode by mode:
#
#     ⟨ψ₀ | e^{-iH_f t} | ψ₀⟩
#       = ∏_k  ⟨n₀(k) | e^{-iε_{J_f}(k) t (n̂_k - 1/2)} | n₀(k)⟩
#       = ∏_k  exp{ -i ε_{J_f}(k) t · (n₀(k) - 1/2) },
#
# where `n₀(k) = Θ(-ε_{J₀}(k))` is the initial-state occupation.  The
# **modulus** is
#
#     |⟨ψ₀ | e^{-iH_f t} | ψ₀⟩| = ∏_k 1 = 1                         (♣)
#
# at every site, hence the Loschmidt rate function
#
#     λ(t) = - lim_{N→∞} (1/N) log |⟨ψ₀ | e^{-iH_f t} | ψ₀⟩|² = 0     (♠)
#
# is identically **zero** whenever the Fermi sea of the initial and
# the final Hamiltonian agree (`sgn J₀ = sgn J_f`).
#
# ── Sign-flip quench (orthogonality) ─────────────────────────────────────────
#
# If `sgn J₀ ≠ sgn J_f` the two Fermi seas are exact complements
# (the modes occupied at J₀ are precisely the ones empty at J_f, and
# vice-versa).  In the thermodynamic limit any two such Slater
# determinants are exactly orthogonal in the strong sense:
#
#     |⟨GS(J₀) | GS(J_f)⟩| = 0,
#
# hence `λ(t) = +∞` for all `t ≥ 0`.  This is the
# Anderson-orthogonality limit — well-defined as a divergent rate but
# numerically degenerate.  We expose it with a warning rather than a
# finite answer.
#
# ── Why the Loschmidt rate of an XX → XX quench is degenerate ────────────────
#
# Reference: Calabrese, Essler & Fagotti, J. Stat. Mech. (2012) P07016
# treat XX-quench Loschmidt **from a Néel / dimer initial state** —
# states that are *not* Gaussian in the same fermion basis and
# therefore induce a non-trivial single-particle Bogoliubov rotation
# at the quench instant.  The XXZ1D model in QAtlas has no
# magnetic-field, dimerisation, or Néel-state machinery yet, so the
# only XX → XX quench expressible in the current model class is
# `|GS(J₀)⟩ → e^{-iH_XX(J_f)t} |GS(J₀)⟩`, which is the degenerate
# `(♣) / (♠)` case derived above.
#
# Phase-1 deliverable: closed-form `λ(t) = 0` (same-sign J quench),
# and `+∞` (sign-flip) with rigorous documentation.  Phase 2
# (deferred) will add either (a) a magnetic-field generalisation of
# `XXZ1D`, or (b) a separate `XYModel` carrying the pairing γ, at
# which point the Loschmidt rate becomes the textbook
# Calabrese-Essler-Fagotti integral.
#
# References
#   - P. Calabrese, F.H.L. Essler, M. Fagotti, J. Stat. Mech. (2012)
#     P07016.
#   - M. Heyl, A. Polkovnikov, S. Kehrein, Phys. Rev. Lett. 110,
#     135704 (2013) — Loschmidt rate / dynamical phase transitions.
#   - F.H.L. Essler, M. Fagotti, J. Stat. Mech. (2016) 064002 — review
#     of quench in integrable systems.
# ─────────────────────────────────────────────────────────────────────────────

# ─── XX-quench guard ────────────────────────────────────────────────────────
#
# The same `_xx_is_free_fermion` predicate as PR #157 (`XXZ_xx_infinite.jl`)
# is *not* assumed to be in scope yet — phase-1 lands on `main`, phase-1.5
# lands `XXZ_xx_infinite.jl`, and the include order in src/QAtlas.jl
# is independent.  Define a phase-148-local copy and gate every
# Δ-sensitive entry point through it.
@inline _xx_quench_is_free_fermion(model::XXZ1D) =
    isapprox(model.Δ, 0.0; atol=1e-12)

# Throw a precise DomainError so callers (and `iszero(Δ)` filter
# tests) can match on the type.
function _xx_quench_assert_free_fermion(model_f::XXZ1D, model_0::XXZ1D)
    if !_xx_quench_is_free_fermion(model_f)
        throw(
            DomainError(
                model_f.Δ,
                "XXZ1D LoschmidtEcho is implemented only at the XX point " *
                "Δ = 0 in this release (issue #148 phase 1).  The general-Δ " *
                "thermal Bethe-ansatz route (issue #108) is a separate workstream.",
            ),
        )
    end
    if !_xx_quench_is_free_fermion(model_0)
        throw(
            DomainError(
                model_0.Δ,
                "XXZ1D LoschmidtEcho: the *initial* model must also be at " *
                "Δ = 0 in this release; got initial.Δ = $(model_0.Δ).",
            ),
        )
    end
    return nothing
end

# ─── Same-Fermi-sea predicate ────────────────────────────────────────────────
#
# In the spin convention `ε_J(k) = J cos(k)` (XXZ.jl §★), so
# `n_k(J) = Θ(-J cos k)`.  Since `cos k > 0` on (-π/2, π/2) and
# `cos k < 0` on (π/2, π) ∪ (-π, -π/2), the Fermi sea
#
#     {k : J cos k < 0}
#
# depends on `sgn J` only.  Two J's with the same sign therefore share
# the Fermi sea; two J's with opposite signs have disjoint
# (complementary) Fermi seas; J = 0 is the trivial flat-band limit
# (every k is degenerate at zero energy).
#
# The function returns
#   :same      if sgn J₀ == sgn J_f != 0          (degenerate quench)
#   :flipped   if sgn J₀ == -sgn J_f != 0         (orthogonality)
#   :flat_initial   if J₀ == 0 (initial Hamiltonian is the flat band)
#   :flat_final     if J_f == 0 (final Hamiltonian is the flat band)
function _xx_quench_fermi_sea_class(J0::Real, Jf::Real)::Symbol
    if iszero(J0) && iszero(Jf)
        return :flat_both
    elseif iszero(J0)
        return :flat_initial
    elseif iszero(Jf)
        return :flat_final
    elseif sign(J0) == sign(Jf)
        return :same
    else
        return :flipped
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Public fetch dispatch for `LoschmidtEcho{:rate}` at Infinite().
#
#   fetch(model_f::XXZ1D, ::LoschmidtEcho{:rate}, ::Infinite;
#         initial::XXZ1D, t::Real) -> Float64
#
# Δ = 0 only.  Returns
#
#   λ(t) = 0                  for sgn J₀ == sgn J_f != 0
#                              (Fermi seas coincide; the amplitude is
#                              a pure phase, |L(t)| = 1)
#   λ(t) = Inf                for sgn J₀ != sgn J_f, both nonzero
#                              (Anderson orthogonality; emits @warn)
#   λ(t) = 0                  for the J₀ = J_f = 0 flat-band fixed point
#                              (the state is annihilated by H_f, no
#                              dynamics whatsoever)
#   λ(t) = NaN + @warn        for the mixed flat-band cases (J₀=0, J_f≠0
#                              or vice-versa); the GS of a flat band is
#                              degenerate so a definite Loschmidt rate
#                              requires an initial-state choice the
#                              issue does not provide.
#
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(model_f::XXZ1D, ::LoschmidtEcho{:rate}, ::Infinite;
          initial::XXZ1D, t::Real) -> Float64

Loschmidt rate function

    λ(t) = - lim_{N→∞} (1/N) log |⟨ψ₀ | e^{-iH_f t} | ψ₀⟩|²

for the **XX → XX** quench `H_XX(J₀) → H_XX(J_f)` of the infinite
Δ = 0 chain (Calabrese-Essler-Fagotti, J. Stat. Mech. (2012) P07016,
specialised to the Gaussian-state-to-Gaussian-state, no-pairing
sub-case).

Because `H_XX(J)` Jordan–Wigner-transforms to a tight-binding fermion
chain *without* pairing, both `H_XX(J₀)` and `H_XX(J_f)` are
diagonalised in the **same** plane-wave basis `c(k)`.  The Loschmidt
amplitude therefore factorises into single-mode phases,

    ⟨ψ₀ | e^{-iH_f t} | ψ₀⟩
      = ∏_k exp{ -i ε_{J_f}(k) t · (n_k(J₀) - 1/2) },

whose modulus is **identically 1**.  Hence

    λ(t) = 0      whenever sgn J₀ = sgn J_f                 (same Fermi sea)

and a divergent λ(t) = +∞ (Anderson orthogonality of the two Fermi
seas) when `sgn J₀ ≠ sgn J_f`.

# Arguments

- `model_f::XXZ1D` — final Hamiltonian.  Must have `Δ = 0` (otherwise
  `DomainError`); a follow-up issue (#108 / #143) covers `Δ ≠ 0`.
- `initial::XXZ1D` — initial Hamiltonian whose ground state is the
  pre-quench state `|ψ₀⟩`.  Must also have `Δ = 0`.
- `t::Real` — real evolution time.

# Returns

A `Float64` with `λ(t) ≥ 0`; degenerate cases as documented above.

# Examples

```julia-repl
julia> m = XXZ1D(; J=1.0, Δ=0.0);

julia> fetch(m, LoschmidtEcho(; mode=:rate), Infinite(); initial=m, t=1.0)
0.0

julia> fetch(m, LoschmidtEcho(; mode=:rate), Infinite();
             initial=XXZ1D(; J=0.5, Δ=0.0), t=1.0)
0.0
```

# References

- P. Calabrese, F.H.L. Essler, M. Fagotti, *J. Stat. Mech.* (2012)
  P07016 — XX-quench dynamics, Loschmidt amplitude.
- M. Heyl, A. Polkovnikov, S. Kehrein, *Phys. Rev. Lett.* 110, 135704
  (2013) — definition of the dynamical Loschmidt rate λ(t).
- F.H.L. Essler, M. Fagotti, *J. Stat. Mech.* (2016) 064002.
"""
function fetch(
    model_f::XXZ1D,
    ::LoschmidtEcho{:rate},
    ::Infinite;
    initial::XXZ1D,
    t::Real,
    kwargs...,
)
    _xx_quench_assert_free_fermion(model_f, initial)

    cls = _xx_quench_fermi_sea_class(initial.J, model_f.J)
    if cls === :same
        # ε_{J_f}(k) is real, the Fermi sea coincides with the initial
        # one, the amplitude is a pure phase, and the rate vanishes
        # exactly.  See `(♣) / (♠)` in the file header.
        return 0.0
    elseif cls === :flipped
        # Anderson orthogonality of complementary Fermi seas.  λ = +∞.
        @warn (
            "XXZ1D LoschmidtEcho{:rate}: sgn J₀ ≠ sgn J_f at the XX " *
            "point makes the initial and final Fermi seas complementary " *
            "(Anderson orthogonality).  The Loschmidt rate is +∞ for " *
            "every t ≥ 0; returning Inf."
        ) J_initial = initial.J J_final = model_f.J
        return Inf
    elseif cls === :flat_both
        # Both Hamiltonians are zero — no dynamics at all.
        return 0.0
    else  # :flat_initial or :flat_final
        # The flat-band ground state is highly degenerate; the Loschmidt
        # rate then depends on which member of the degenerate manifold
        # is taken as |ψ₀⟩, which the issue scope does not specify.
        @warn (
            "XXZ1D LoschmidtEcho{:rate}: one of {J_initial, J_final} is " *
            "zero (flat band).  The flat-band ground state is " *
            "exponentially degenerate and the rate is not single-valued; " *
            "returning NaN.  Provide both J_initial ≠ 0 and J_final ≠ 0 " *
            "for a well-defined answer."
        ) J_initial = initial.J J_final = model_f.J
        return NaN
    end
end
