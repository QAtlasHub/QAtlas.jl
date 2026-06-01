# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE — Juttner-Klumper-Suzuki QTM NLIE module (#523)
#
# Phase-2B foundation for the full thermodynamics of the 1D Hubbard model
# at arbitrary (T, U, mu, h). Stage A (this file) ships:
#
#   (i)  the atomic-limit closed-form free energy (t = 0), exact for any
#        (T, U, mu, h)
#   (ii) a typed Grid skeleton (`JKSGrid`) that Stage B will fill with
#        full NLIE contour + auxiliary-function storage
#   (iii) stub signatures for the kernels and driving functions Stage B
#        will implement
#
# No NLIE solver is wired up in Stage A. Production fetch dispatches for
# `Hubbard1D` still go through the regime-based stopgap in
# `Hubbard1D_thermal_stopgap.jl` (PR #530). Stage B will route them
# through the full Picard NLIE solver and remove the NaN intermediate
# regime.
#
# Math reference: G. Juttner, A. Klumper, J. Suzuki, "The Hubbard chain
# at finite temperatures: ab initio calculations of Tomonaga-Luttinger
# liquid properties", Nucl. Phys. B 522, 471 (1998), arXiv:cond-mat/
# 9711310. Numbering below cites that paper's equations.
# ─────────────────────────────────────────────────────────────────────────────

module Hubbard1DJKSNLIE

export atomic_free_energy, atomic_free_energy_half_filling

# ═══════════════════════════════════════════════════════════════════════════════
# Atomic limit (t = 0) — closed form
# ═══════════════════════════════════════════════════════════════════════════════
#
# At t = 0 the Hubbard sites decouple. The grand-canonical per-site
# partition function for U, mu, H (H = magnetic field) is
#
#     Z_site = sum over {empty, up, down, doubly} of e^{-beta(E - mu N - h M)}
#            = 1 + e^{beta(mu + h/2)} + e^{beta(mu - h/2)} + e^{beta(2 mu - U)}
#
# yielding f_atomic = -T ln Z_site. At H = 0 this collapses to
#
#     Z_site = 1 + 2 e^{beta mu} + e^{beta(2 mu - U)}
#
# which drives the JKS high-T asymptotic analysis (Sec 7.3). Limits:
#
#   - T -> infinity, any (U, mu, h):  f -> -T ln 4   (all 4 states equally weighted)
#   - T -> 0,        mu = U/2:        f -> -U/2 - T ln 2 (singly-occupied doublet
#                                                       residual spin entropy)

"""
    atomic_free_energy(beta, U, mu; h=0.0) -> Float64

Atomic-limit (t = 0) Helmholtz free-energy density of the Hubbard model
in the grand canonical ensemble. Returns

    f = -T ln[1 + e^{beta(mu + h/2)} + e^{beta(mu - h/2)} + e^{beta(2 mu - U)}].

Exact for any (T, U, mu, h) at t = 0 and reproduces JKS Sec 7.3 high-T
limit when used as the leading-order term in beta * t.

Throws `DomainError` for `beta <= 0`.
"""
function atomic_free_energy(beta::Real, U::Real, mu::Real; h::Real=0.0)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    Z = 1 + exp(beta * (mu + h/2)) + exp(beta * (mu - h/2)) + exp(beta * (2 * mu - U))
    return -log(Z) / beta
end

"""
    atomic_free_energy_half_filling(beta, U) -> Float64

Specialization of `atomic_free_energy` at half filling (`mu = U/2`,
`h = 0`):

    f = -T ln[2(1 + e^{beta U / 2})] = -T ln 2 - T ln(1 + e^{beta U / 2}).

Reduces to `-T ln 4` at `U = 0` and `-U/2 - T ln 2` at low T (singly-
occupied doublet residual entropy).
"""
function atomic_free_energy_half_filling(beta::Real, U::Real)
    return atomic_free_energy(beta, U, U/2)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Stage-B placeholders — NLIE grid + kernels + driving functions
# ═══════════════════════════════════════════════════════════════════════════════
#
# The structs and signatures below establish the API that Stage B's
# full Picard solver will fill. They error in Stage A to make accidental
# calls during this phase loud rather than silent.

"""
    JKSGrid(N, gamma)

NLIE discretization grid for the JKS QTM NLIE (Stage B). At full
implementation will hold the contour points, kernel cache, and
auxiliary-function arrays for `b, b_bar, c, c_bar`. Stage A reserves
the type for forward compatibility only.

# Fields

- `N::Int` — number of contour points (default 128 in Stage B)
- `gamma::Float64` — contour-shift parameter `gamma in (0, pi)` (eq 51,
  default `2*pi/3`)
"""
struct JKSGrid
    N::Int
    gamma::Float64
    function JKSGrid(N::Int, gamma::Real)
        N > 0 || throw(DomainError(N, "N must be > 0"))
        0 < gamma < pi || throw(DomainError(gamma, "gamma must be in (0, pi)"))
        return new(N, Float64(gamma))
    end
end

"""
    jks_kernel_K_n(s, n, gamma) -> ComplexF64

JKS kernel `K_n(s) = gamma / (pi * s * (s + 2 n i gamma))` from eq (38),
used in the NLIE convolutions of Sec 5. Stage A: not yet implemented;
calling it raises `ErrorException` so accidental Stage-A use is loud.

Stage B will provide the full evaluator and a cached Toeplitz/full
convolution operator.
"""
function jks_kernel_K_n(s::Number, n::Integer, gamma::Real)
    error("Hubbard1DJKSNLIE.jks_kernel_K_n is not yet implemented (Stage B / #523).")
end

"""
    jks_driving_b(x; h=0.0) -> Float64

JKS spin-channel driving term `d_b(x) = -h` from eq (54). Constant in
`x` because the field couples uniformly to the magnetization.

This one is fully implemented (it is trivial); the corresponding c /
c_bar driving terms still need the Trotter-limit function
`log lambda_0(x)` that Stage B will derive from Sec 4.
"""
jks_driving_b(x::Real; h::Real=0.0) = -float(h)

"""
    jks_driving_c(x; U, mu, h=0.0) -> Float64

JKS charge-channel driving term `d_c(x) = -U/2 + (mu + h/2) + log lambda_0(x)`
from eq (54). Stage A returns the constant part only; the Trotter-limit
`log lambda_0(x)` term is deferred to Stage B. Calling Stage A is loud:
we error so a Stage-A caller does not silently get a wrong value.
"""
function jks_driving_c(x::Real; U::Real, mu::Real, h::Real=0.0)
    error(
        "Hubbard1DJKSNLIE.jks_driving_c needs Trotter-limit log lambda_0; deferred to Stage B (#523).",
    )
end

"""
    jks_driving_cbar(x; U, mu, h=0.0) -> Float64

Conjugate of `jks_driving_c`: `d_c_bar(x) = -U/2 - (mu + h/2) - log lambda_0(x)`
from eq (55). Stage A error; full implementation in Stage B.
"""
function jks_driving_cbar(x::Real; U::Real, mu::Real, h::Real=0.0)
    error(
        "Hubbard1DJKSNLIE.jks_driving_cbar needs Trotter-limit log lambda_0; deferred to Stage B (#523).",
    )
end

end  # module Hubbard1DJKSNLIE
