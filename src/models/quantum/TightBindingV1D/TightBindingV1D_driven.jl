# ─────────────────────────────────────────────────────────────────────────────
# TightBindingV1D_driven.jl — exact nonlinear response of the ac-driven
# free-fermion (V = 0) tight-binding chain.
#
# A spatially-uniform monochromatic electric field E(τ) = E₀ cos(ωτ) couples to
# the chain through the Peierls substitution: the hopping acquires a
# time-dependent phase t → t e^{±iA(τ)} with the vector potential (e = ℏ = a = 1)
#
#     A(τ) = -∫ E dτ = -(E₀/ω) sin(ωτ) ≡ -K sin(ωτ),     K ≡ E₀/ω  (dimensionless).
#
# The intraband group-velocity current of a Bloch mode at crystal momentum k is
# ε'(k) evaluated at the shifted momentum k - A(τ):
#
#     j(k, τ) = 2t sin(k - A(τ)) = 2t sin(k + K sin ωτ).
#
# Jacobi–Anger (Abramowitz–Stegun 9.1.42-43 / DLMF 10.12) expands this into an
# EXACT harmonic series whose amplitudes are Bessel functions of the first kind:
#
#     j(k, τ) = 2t J₀(K) sin k                                   (dc / 0-th)
#             + 4t sin k Σ_{m≥1} J_{2m}(K)   cos(2mωτ)           (even harmonics)
#             + 4t cos k Σ_{m≥0} J_{2m+1}(K) sin((2m+1)ωτ)       (odd  harmonics).
#
# Two exact, nonperturbative ("all orders in the field") statements follow, both
# verified against an independent RK4 real-time simulation in
# test/identities/test_driven_free_fermion_nonlinear_response.jl:
#
#   • DYNAMIC LOCALIZATION (Dunlap–Kenkre 1986; Holthaus–Hone 1996).  The dc /
#     cycle-averaged current — hence the coherent transport — is renormalized by
#     the Bessel factor t_eff = t J₀(K).  The band collapses (t_eff → 0) at the
#     zeros of J₀ (first at K = 2.404826…): a static tilt no longer drives a
#     current even though E₀ ≠ 0.  Exposed as `fetch(_, DynamicLocalization, _)`.
#
#   • HARMONIC SPECTRUM.  The n-th harmonic of the current is exactly ∝ Jₙ(K)
#     (`driven_band_harmonic_weights`).  For small K, Jₙ(K) ≈ (K/2)ⁿ/n!, so the
#     n-th harmonic is the order-n (χ⁽ⁿ⁾) nonlinear response and n = 1 recovers
#     the linear conductivity.
#
# NOTE (V = 0 only).  Peierls dressing keeps the chain quadratic only at the
# free-fermion point; V ≠ 0 (Jordan-Wigner-equivalent to the interacting XXZ
# chain) is deferred to Phase 2 and raises `DomainError`, matching the rest of
# this model file.
#
# References:
#   - D. H. Dunlap, V. M. Kenkre, "Dynamic localization of a charged particle
#     moving under the influence of an electric field", Phys. Rev. B 34, 3625
#     (1986).
#   - M. Holthaus, D. W. Hone, "Localization effects in ac-driven tight-binding
#     lattices", [HolthausHone1996](@cite).
# ─────────────────────────────────────────────────────────────────────────────

using SpecialFunctions: besselj, besselj0

"""
    fetch(m::TightBindingV1D, ::DynamicLocalization, ::Infinite;
          drive, t=m.t, V=m.V, kwargs...) -> Float64

Renormalized hopping `t_eff = t · J₀(K)` of the ac-driven free-fermion (V = 0)
tight-binding chain, where `drive = K = E₀/ω` is the dimensionless field
amplitude / frequency (Peierls coupling, e = ℏ = a = 1).

`t_eff` sets the collapsed bandwidth `4 t_eff`; dynamic localization is the
vanishing `t_eff = 0` at the zeros of `J₀` (first at `K ≈ 2.404826`), where a
static tilt drives no current despite `E₀ ≠ 0` (Dunlap–Kenkre 1986;
Holthaus–Hone 1996).  The full harmonic spectrum of the current is
`driven_band_harmonic_weights`.

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::DynamicLocalization,
    ::Infinite;
    drive::Real,
    t::Real=m.t,
    V::Real=m.V,
    kwargs...,
)
    t > 0 || throw(
        DomainError(t, "TightBindingV1D DynamicLocalization requires t > 0; got t = $t."),
    )
    if !iszero(V)
        throw(
            DomainError(
                V,
                "TightBindingV1D DynamicLocalization: V ≠ 0 (JW-equivalent to interacting " *
                "XXZ, Yang-Yang 1966) deferred to Phase 2. Got V = $V.",
            ),
        )
    end
    return t * besselj0(drive)
end

"""
    driven_band_harmonic_weights(drive::Real; nmax::Integer=6) -> Vector{Float64}

Exact Bessel weights `[J₀(K), J₁(K), …, J_nmax(K)]` of the current harmonics of
a single Bloch mode of the ac-driven free-fermion tight-binding chain, with
`drive = K = E₀/ω`.

Entry `n+1` is the Bessel weight of the n-th harmonic `n ω` in the Jacobi–Anger
expansion of `j(k, τ) = 2t sin(k + K sin ωτ)`.  The physical n-th-harmonic
amplitude of the current is `2t·|J₀(K) sin k|` for `n = 0` and
`4t·|Jₙ(K)|·(|sin k| if n even else |cos k|)` for `n ≥ 1`.  For small `K`,
`Jₙ(K) ≈ (K/2)ⁿ/n!`, exhibiting the n-th harmonic as the order-n (χ⁽ⁿ⁾)
nonlinear response.
"""
function driven_band_harmonic_weights(drive::Real; nmax::Integer=6)
    nmax >= 0 || throw(
        DomainError(nmax, "driven_band_harmonic_weights requires nmax ≥ 0; got $nmax.")
    )
    return Float64[besselj(n, drive) for n in 0:nmax]
end

"""
    fetch(m::TightBindingV1D, ::HighHarmonicAmplitude, ::Infinite;
          drive, harmonic, t=m.t, V=m.V, kwargs...) -> Float64

Peak amplitude of the `harmonic`-th harmonic (`n ω`) of the intraband current of
the ac-driven free-fermion (V = 0) chain — the exact, all-orders-in-field
higher-order response.  With `drive = K = E₀/ω`, maximizing over crystal
momentum,

    A₀(K) = 2t |J₀(K)|,      Aₙ(K) = 4t |Jₙ(K)|   (n ≥ 1).

`harmonic = 1` is the linear response; `harmonic ≥ 2` the nonlinear higher
harmonics.  The leading small-field (perturbative χ⁽ⁿ⁾) coefficient is
`nonlinear_susceptibility(; order=harmonic)`.

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::HighHarmonicAmplitude,
    ::Infinite;
    drive::Real,
    harmonic::Integer,
    t::Real=m.t,
    V::Real=m.V,
    kwargs...,
)
    t > 0 || throw(
        DomainError(t, "TightBindingV1D HighHarmonicAmplitude requires t > 0; got t = $t."),
    )
    harmonic >= 0 || throw(
        DomainError(
            harmonic,
            "TightBindingV1D HighHarmonicAmplitude requires harmonic ≥ 0; got $harmonic.",
        ),
    )
    if !iszero(V)
        throw(
            DomainError(
                V,
                "TightBindingV1D HighHarmonicAmplitude: V ≠ 0 (JW-equivalent to interacting " *
                "XXZ, Yang-Yang 1966) deferred to Phase 2. Got V = $V.",
            ),
        )
    end
    jn = driven_band_harmonic_weights(drive; nmax=harmonic)[harmonic + 1]  # Jₙ(K)
    return (harmonic == 0 ? 2 : 4) * t * abs(jn)
end

"""
    nonlinear_susceptibility(; order::Integer, omega::Real=1.0, t::Real=1.0) -> Float64

Leading-order (perturbative χ⁽ⁿ⁾) nonlinear conductivity of the ac-driven
free-fermion band: the coefficient of `E₀^order` in the peak `order`-th harmonic
current, obtained from the small-drive limit `Jₙ(K) ≈ (K/2)ⁿ/n!` (`K = E₀/ω`):

    χ⁽ⁿ⁾ = 4t / (n! (2ω)ⁿ),     n = order ≥ 1.

`order = 1` is the linear (Drude) conductivity `2t/ω`; `order = 2` the
rectification / second-harmonic coefficient `t/(2ω²)`; and so on.  The exact,
all-orders amplitude is [`HighHarmonicAmplitude`](@ref) = `4t|Jₙ(E₀/ω)|`, of
which this is the `K → 0` leading term.
"""
function nonlinear_susceptibility(; order::Integer, omega::Real=1.0, t::Real=1.0)
    order >= 1 || throw(
        DomainError(order, "nonlinear_susceptibility requires order ≥ 1; got $order.")
    )
    omega > 0 ||
        throw(DomainError(omega, "nonlinear_susceptibility requires ω > 0; got $omega."))
    # BigInt/BigFloat denominator so large `order` cannot silently overflow Int64.
    return Float64(4 * t / (factorial(big(order)) * big(2 * omega)^order))
end
