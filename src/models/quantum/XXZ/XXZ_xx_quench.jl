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
# where `n₀(k) = Θ(-ε_{J₀}(k))` is the initial-state occupation.  Every
# factor in the product is a **pure phase** (modulus 1), so the
# **modulus** of the full amplitude is
#
#     |⟨ψ₀ | e^{-iH_f t} | ψ₀⟩| = ∏_k 1 = 1                         (♣)
#
# **independently of the signs of J₀ and J_f**, hence the Loschmidt rate
# function
#
#     λ(t) = - lim_{N→∞} (1/N) log |⟨ψ₀ | e^{-iH_f t} | ψ₀⟩|² = 0     (♠)
#
# is identically **zero** for every (J₀, J_f), including the sign-flip
# case `sgn J₀ ≠ sgn J_f`.
#
# ── Why the sign-flip case is *not* Anderson orthogonality ───────────────────
#
# It is true that `|⟨GS(J₀) | GS(J_f)⟩| = 0` when the two Fermi seas are
# complementary (sign-flip): the two Slater determinants have disjoint
# occupied-mode sets and are exactly orthogonal.  But that is the
# **static overlap of two ground states**, which is *not* the Loschmidt
# amplitude.  The Loschmidt amplitude is the autocorrelation
# `⟨ψ₀ | e^{-iH_f t} | ψ₀⟩`, and for the XX → XX quench
# `|ψ₀⟩ = |GS(J₀)⟩` is also a number eigenstate of `H_f`
# (since `H_f = Σ_k ε_{J_f}(k) ĉ_k^† ĉ_k` is diagonal in the same
# plane-wave basis), so
#
#     e^{-iH_f t} |ψ₀⟩ = e^{-iE_f^{(J₀)} t} |ψ₀⟩
#
# with `E_f^{(J₀)} = Σ_{k∈Fermi sea(J₀)} ε_{J_f}(k)`.  At sign-flip the
# state |ψ₀⟩ is the *highest-energy* eigenstate of `H_f` rather than the
# ground state of `H_f`, but it is still an exact eigenstate, so the
# autocorrelation has modulus 1 and λ(t) = 0.
#
# Anderson orthogonality `⟨GS(J₀)|GS(J_f)⟩ = 0` is therefore irrelevant
# to this quench protocol.  (It would matter if the quench observable
# were `|⟨ψ₀|ψ_f⟩|²` — the *initial-to-final-GS* overlap — but the
# Loschmidt rate is by definition an autocorrelation.)
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
# `|GS(J₀)⟩ → e^{-iH_XX(J_f)t} |GS(J₀)⟩`, which is the trivial
# `λ(t) ≡ 0` case derived above.
#
# Phase-1 deliverable: closed-form `λ(t) ≡ 0` for the GS-to-GS XX → XX
# quench (all sgn-J combinations).  Phase 2 (deferred) will add either
# (a) a magnetic-field generalisation of `XXZ1D`, or (b) a separate
# `XYModel` carrying the pairing γ, at which point the Loschmidt rate
# becomes the textbook Calabrese-Essler-Fagotti integral.
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
@inline _xx_quench_is_free_fermion(model::XXZ1D) = isapprox(model.Δ, 0.0; atol=1e-12)

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

# ─────────────────────────────────────────────────────────────────────────────
# Public fetch dispatch for `LoschmidtEcho{:rate}` at Infinite().
#
#   fetch(model_f::XXZ1D, ::LoschmidtEcho{:rate}, ::Infinite;
#         initial::XXZ1D, t::Real) -> Float64
#
# Δ = 0 only.  Returns `λ(t) ≡ 0` for every (J₀, J_f) — see `(♣) / (♠)`
# in the file header.  The flat-band and sign-flip edge cases all
# reduce to the same `|L(t)| = 1` argument because |ψ₀⟩ is a number
# eigenstate of H_f in the shared plane-wave basis.
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

whose modulus is **identically 1** (each factor is a pure phase).
Equivalently, `|ψ₀⟩ = |GS(J₀)⟩` is a number eigenstate of `H_f` because
both Hamiltonians are diagonal in the same `{n̂_k}`, so
`e^{-iH_f t}|ψ₀⟩ = e^{-iE_f^{(J₀)} t}|ψ₀⟩` is a pure phase.  Hence

    λ(t) ≡ 0          for every (J₀, J_f), including sgn J₀ ≠ sgn J_f.

The static overlap `|⟨GS(J₀)|GS(J_f)⟩|` *does* vanish at sign-flip
(Anderson orthogonality of complementary Fermi seas), but that is a
different quantity from the Loschmidt autocorrelation and does not
enter `λ(t)`.

This trivial result holds because the only XX → XX quench expressible
in the current `XXZ1D` model class is GS-to-GS without any
Bogoliubov-rotation knob.  A non-trivial Loschmidt rate appears in
quenches from non-Gaussian initial states (Néel / dimer; CEF 2012) or
under XY-pairing dynamics, both of which are Phase-2 follow-ups
(see issues #143, #146).

# Arguments

- `model_f::XXZ1D` — final Hamiltonian.  Must have `Δ = 0` (otherwise
  `DomainError`); a follow-up issue (#108 / #143) covers `Δ ≠ 0`.
- `initial::XXZ1D` — initial Hamiltonian whose ground state is the
  pre-quench state `|ψ₀⟩`.  Must also have `Δ = 0`.
- `t::Real` — real evolution time.

# Returns

`0.0` (the trivial Phase-1 value, valid for every `t ≥ 0` and every
sign combination of `(initial.J, model_f.J)`).

# Examples

```julia-repl
julia> m = XXZ1D(; J=1.0, Δ=0.0);

julia> fetch(m, LoschmidtEcho(; mode=:rate), Infinite(); initial=m, t=1.0)
0.0

julia> fetch(m, LoschmidtEcho(; mode=:rate), Infinite();
             initial=XXZ1D(; J=0.5, Δ=0.0), t=1.0)
0.0

julia> fetch(m, LoschmidtEcho(; mode=:rate), Infinite();
             initial=XXZ1D(; J=-1.0, Δ=0.0), t=1.0)   # sign-flip — also 0
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
    model_f::XXZ1D, ::LoschmidtEcho{:rate}, ::Infinite; initial::XXZ1D, t::Real, kwargs...
)
    _xx_quench_assert_free_fermion(model_f, initial)
    # |ψ₀⟩ = |GS(initial.J)⟩ is a number eigenstate of H_XX(model_f.J)
    # in the shared plane-wave basis, so the Loschmidt amplitude is a
    # pure phase, |L(t)| = 1, and λ(t) ≡ 0 for every (J₀, J_f, t).
    # See `(♣) / (♠)` in the file header for the full derivation.
    return 0.0
end

# -----------------------------------------------------------------------------
# Lieb-Robinson velocity + entanglement growth slope at the XX point (Delta = 0)
# -----------------------------------------------------------------------------

"""
    fetch(model::XXZ1D, ::LiebRobinsonVelocity, ::Infinite;
          J = model.J, Delta = model.Δ, kwargs...) -> Float64

Maximum group velocity of the free-fermion XX chain (`Delta = 0`).
The single-particle dispersion is `epsilon(k) = -2 J cos k`, and its
group-velocity maximum is

    v_LR = 2 |J|,

attained at `k = pi / 2`. For `Delta != 0` the LR bound is more
involved (interacting regime) and is deferred — this dispatch throws
`DomainError`.
"""
function fetch(
    model::XXZ1D,
    ::LiebRobinsonVelocity,
    ::Infinite;
    J::Real=model.J,
    Delta::Real=model.Δ,
    kwargs...,
)
    isapprox(Delta, 0.0; atol=1e-12) || throw(
        DomainError(
            Delta,
            "XXZ1D LiebRobinsonVelocity: closed form only at the XX point (Delta = 0); got Delta = $Delta.",
        ),
    )
    return 2 * abs(J)
end

"""
    fetch(model::XXZ1D, ::EntanglementGrowthSlope, ::Infinite;
          beta_eff::Real, kwargs...) -> Float64

Calabrese-Cardy 2005 linear-growth slope of half-system entanglement
after a quench at the XX point of the XXZ chain. The chain is
gapless with central charge `c = 1` and LR velocity `v = 2 |J|`, so

    dS_A / dt = pi c v / (3 beta_eff) = 2 pi |J| / (3 beta_eff).

Delegates to `Universality(:XY)` (c = 1).
"""
function fetch(
    model::XXZ1D, ::EntanglementGrowthSlope, ::Infinite; beta_eff::Real, kwargs...
)
    isapprox(model.Δ, 0.0; atol=1e-12) || throw(
        DomainError(
            model.Δ,
            "XXZ1D EntanglementGrowthSlope: closed form only at the XX point (Delta = 0); got Delta = $(model.Δ).",
        ),
    )
    return fetch(
        Universality(:XY),
        EntanglementGrowthSlope(),
        Infinite();
        v=fetch(model, LiebRobinsonVelocity(), Infinite()),
        beta_eff=beta_eff,
        kwargs...,
    )
end

"""
    fetch(model::XXZ1D, ::EntanglementSaturationDensity, ::Infinite;
          beta_eff::Real, kwargs...) -> Float64

Long-time saturation of post-quench half-system entanglement entropy
at the XX point (`Delta = 0`). The XX chain is gapless with c = 1,
so delegates to `Universality(:XY)` returning `pi / (6 beta_eff)`.

Off-XX (`Delta != 0`) throws `DomainError` -- interacting regime
deferred.
"""
function fetch(
    model::XXZ1D, ::EntanglementSaturationDensity, ::Infinite; beta_eff::Real, kwargs...
)
    isapprox(model.Δ, 0.0; atol=1e-12) || throw(
        DomainError(
            model.Δ,
            "XXZ1D EntanglementSaturationDensity: closed form only at the XX " *
            "point (Delta = 0); got Delta = $(model.Δ).",
        ),
    )
    return fetch(
        Universality(:XY),
        EntanglementSaturationDensity(),
        Infinite();
        beta_eff=beta_eff,
        kwargs...,
    )
end
