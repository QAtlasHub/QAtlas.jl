# models/quantum/AKLT/AKLT1D_thermal.jl — biquadratic-aware high-temperature
# series expansion (HTSE) of the infinite AKLT chain's thermodynamics (#506).
#
# The β = ∞ exact limits live in AKLT1D.jl as the `scheme = :canonical`
# definitions of (AKLT1D, {FreeEnergy, SpecificHeat, ThermalEntropy}, Infinite).
# Here we add the `scheme = :htse` approximation that covers the *finite*-β
# high-temperature regime those canonical rows reject.
#
# Method.  The per-site log-partition function has the linked-cluster
# (cumulant) expansion
#
#     φ(β) ≡ ln Z / N = ln 3 + Σ_{n≥1} aₙ βⁿ ,   aₙ = (−1)ⁿ κₙ Jⁿ / n! ,
#
# where κₙ are the per-bond infinite-temperature cumulants of the AKLT bond
# operator  h_b = S_i·S_{i+1} + (1/3)(S_i·S_{i+1})²  (J = 1).  Each κₙ is the
# thermodynamic-limit increment κₙ = lim_{N→∞}[Cₙ(N) − Cₙ(N−1)] of the
# extensive cumulant Cₙ(N) of H_N, obtained from the exact infinite-T moments
# Tr(Hᵏ)/3ᴺ.  The biquadratic term enters automatically — it is part of h_b —
# which is exactly what the published bilinear-only HTSE machinery
# (Lohmann–Schmidt–Richter 2014, [Lohmann2014]) does not supply.
#
# The four cumulants below are exact rationals (verified by exact-arithmetic ED
# in the test).  Provenance / cross-checks (test/models/quantum/misc/test_aklt_htse.jl):
#   • κ₂ in the *bilinear* limit (drop the biquadratic) equals r²/3 = 4/3 with
#     r = s(s+1) = 2 — the per-bond specific-heat coefficient d₂ of
#     [Lohmann2014] (kagome d₂ = (2/3)r², 2 bonds/site).  This anchors the spin
#     normalisation and the per-site / β conventions against the paper.
#   • the assembled f, c_v, s match an independent Boltzmann-ED of a finite
#     PBC ring to <1% for βJ ≲ 0.35 (the validity window below).
#
# κₙ for n = 1..4 (per bond, J = 1):  4/9, 80/81, −160/729, −6800/2187.
const _AKLT_HTSE_KAPPA = (4 / 9, 80 / 81, -160 / 729, -6800 / 2187)

# Order of the truncation (highest power of β kept in φ).
const _AKLT_HTSE_ORDER = 4

"""
    _aklt_htse_phi_coeffs(J) -> NTuple{4,Float64}

Coefficients `(a₁, a₂, a₃, a₄)` of `φ(β) = ln Z/N = ln 3 + Σ aₙ βⁿ` for the
AKLT chain with coupling `J`, from the per-bond cumulants
[`_AKLT_HTSE_KAPPA`] via `aₙ = (−1)ⁿ κₙ Jⁿ / n!`.
"""
@inline function _aklt_htse_phi_coeffs(J::Real)
    κ = _AKLT_HTSE_KAPPA
    return (-κ[1] * J, κ[2] * J^2 / 2, -κ[3] * J^3 / 6, κ[4] * J^4 / 24)
end

function _aklt_htse_require_beta(beta)
    return (isfinite(beta) && beta > 0) || throw(
        DomainError(
            beta,
            "AKLT1D scheme=:htse is a finite-β high-temperature expansion; needs 0 < β < ∞ (got β = $beta). " *
            "Use scheme=:canonical for the exact β = ∞ limit.",
        ),
    )
end

# Generic fallback: unknown scheme (or a quantity without an :htse definition,
# e.g. SusceptibilityZZ — its HTSE needs separate magnetisation cumulants).
function _aklt_thermo_infinite_scheme(::AKLT1D, q, ::Val{S}; beta) where {S}
    throw(
        ArgumentError(
            "AKLT1D $(typeof(q)) Infinite: no scheme :$(S) (only :canonical + :htse)"
        ),
    )
end

"""
    _aklt_thermo_infinite_scheme(m, ::FreeEnergy, ::Val{:htse}; beta) -> Float64

Per-site Helmholtz free energy `f = −φ(β)/β` of the infinite AKLT chain from
the order-4 biquadratic-aware HTSE.  Valid for `βJ ≲ 0.4` (high temperature).
"""
function _aklt_thermo_infinite_scheme(m::AKLT1D, ::FreeEnergy, ::Val{:htse}; beta)
    _aklt_htse_require_beta(beta)
    a = _aklt_htse_phi_coeffs(m.J)
    φ = log(3) + sum(a[n] * beta^n for n in 1:_AKLT_HTSE_ORDER)
    return -φ / beta
end

"""
    _aklt_thermo_infinite_scheme(m, ::SpecificHeat, ::Val{:htse}; beta) -> Float64

Per-site heat capacity `c_v = Σ n(n−1) aₙ βⁿ` of the infinite AKLT chain from
the order-4 biquadratic-aware HTSE.  Valid for `βJ ≲ 0.4` (high temperature).
"""
function _aklt_thermo_infinite_scheme(m::AKLT1D, ::SpecificHeat, ::Val{:htse}; beta)
    _aklt_htse_require_beta(beta)
    a = _aklt_htse_phi_coeffs(m.J)
    return sum(n * (n - 1) * a[n] * beta^n for n in 2:_AKLT_HTSE_ORDER)
end

"""
    _aklt_thermo_infinite_scheme(m, ::ThermalEntropy, ::Val{:htse}; beta) -> Float64

Per-site Gibbs entropy `s = φ(β) + β ε(β)` (with `ε = −∂φ/∂β`) of the infinite
AKLT chain from the order-4 biquadratic-aware HTSE.  `s → ln 3` as `β → 0`.
Valid for `βJ ≲ 0.4` (high temperature).
"""
function _aklt_thermo_infinite_scheme(m::AKLT1D, ::ThermalEntropy, ::Val{:htse}; beta)
    _aklt_htse_require_beta(beta)
    a = _aklt_htse_phi_coeffs(m.J)
    φ = log(3) + sum(a[n] * beta^n for n in 1:_AKLT_HTSE_ORDER)
    ε = -sum(n * a[n] * beta^(n - 1) for n in 1:_AKLT_HTSE_ORDER)
    return φ + beta * ε
end
