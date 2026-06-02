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

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage B — kernels + discrete convolution scaffold
#
# What this file adds (appended to Hubbard1D_jks_nlie.jl Stage A):
#
#   (1) Concrete `jks_kernel_K_n_concrete(s, n, gamma) -> ComplexF64`
#       evaluator. The Stage-A stub `jks_kernel_K_n` is kept for
#       backward compatibility but is now redirected to the concrete
#       version.
#   (2) `JKSContourGrid(N, gamma; x_max)` struct holding the discretized
#       contour points x_j.
#   (3) `build_kernel_matrix(grid, n) -> Matrix{ComplexF64}` Toeplitz-
#       style discrete convolution operator.
#   (4) `apply_kernel(K_mat, f) -> Vector{ComplexF64}` discrete convolution.
#
# Stage C will use these to write the Picard NLIE residual and solver,
# and to switch the Hubbard1D fetch dispatch to the JKS path.
#
# Reference: JKS 1998 eq (38) for the kernels, eq (50) for the
# discrete-convolution convention.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Kernel — concrete implementation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks_kernel_K_n_concrete(s, n, gamma) -> ComplexF64

JKS kernel `K_n(s) = gamma / (pi * s * (s + 2 n i gamma))` from eq (38).

Vanishes at infinity (`K_n(s) ~ gamma/(pi s^2)` for `|s| -> infty`) and
has simple poles at `s = 0` and `s = -2 n i gamma`.

For `s == 0` returns `Inf + 0im` since the pole at the origin is the
singular feature that the contour deformation in eq (47) exploits;
numerical callers should evaluate at `s + i*eps` (small) instead.
"""
function jks_kernel_K_n_concrete(s::Number, n::Integer, gamma::Real)
    n > 0 || throw(DomainError(n, "kernel index n must be > 0"))
    0 < gamma < pi || throw(DomainError(gamma, "gamma must be in (0, pi)"))
    if s == zero(s)
        return ComplexF64(Inf, 0.0)
    end
    return ComplexF64(gamma / (pi * s * (s + 2 * n * im * gamma)))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Contour grid + kernel-matrix cache
# ═══════════════════════════════════════════════════════════════════════════════

"""
    JKSContourGrid(N, gamma; x_max=10.0)

Discretization of the real-axis contour used by the JKS NLIE. Holds:

- `N::Int` — number of contour points
- `gamma::Float64` — kernel shift parameter, gamma in (0, pi)
- `x_max::Float64` — half-width of the discretized window [-x_max, x_max]
- `x::Vector{Float64}` — N evenly-spaced contour points
- `dx::Float64` — grid spacing (uniform)

The grid is symmetric about 0 to expose the K_n / K_n_bar conjugation
symmetry numerically.
"""
struct JKSContourGrid
    N::Int
    gamma::Float64
    x_max::Float64
    x::Vector{Float64}
    dx::Float64
    function JKSContourGrid(N::Int, gamma::Real; x_max::Real=10.0)
        N > 1 || throw(DomainError(N, "JKSContourGrid requires N > 1"))
        0 < gamma < pi || throw(DomainError(gamma, "gamma must be in (0, pi)"))
        x_max > 0 || throw(DomainError(x_max, "x_max must be > 0"))
        x = collect(range(-x_max, x_max; length=N))
        dx = x[2] - x[1]
        return new(N, Float64(gamma), Float64(x_max), x, dx)
    end
end

"""
    build_kernel_matrix(grid, n; eps=1e-10) -> Matrix{ComplexF64}

Build the discrete convolution matrix
`K[j, k] = jks_kernel_K_n_concrete(x_j - x_k + i*eps, n, gamma) * dx`
so that `(K * f)[j] ≈ integral K_n(x_j - y) f(y) dy` for any vector
`f` sampled on the contour grid. The `eps` (default `1e-10`) shifts
off the singular `x = 0` diagonal — needed because the kernel has a
pole at the origin.
"""
function build_kernel_matrix(grid::JKSContourGrid, n::Integer; eps::Real=1e-10)
    n > 0 || throw(DomainError(n, "kernel index n must be > 0"))
    eps > 0 || throw(DomainError(eps, "eps regularization must be > 0"))
    N = grid.N
    K = zeros(ComplexF64, N, N)
    for j in 1:N, k in 1:N
        K[j, k] =
            jks_kernel_K_n_concrete(grid.x[j] - grid.x[k] + im * eps, n, grid.gamma) *
            grid.dx
    end
    return K
end

"""
    apply_kernel(K, f) -> Vector{ComplexF64}

Discrete convolution `(K * f)[j] = sum_k K[j, k] * f[k]`. Returns a
`ComplexF64` vector. Thin wrapper for readability; matrix-vector
multiply is the same operation.
"""
function apply_kernel(K::AbstractMatrix, f::AbstractVector)
    size(K, 2) == length(f) ||
        throw(DimensionMismatch("K is $(size(K)) but f has length $(length(f))"))
    return K * f
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.1 — auxiliary-function container
#
# Holds the 4 complex-valued auxiliary functions (b, b_bar, c, c_bar) of
# JKS eq (27) on the discretised contour. Provides safe initialisation
# from the atomic-limit fugacities (Stage A) — a reasonable Picard
# starting point for any (beta, U, mu, h), but NOT claimed to be an
# NLIE solution outside the atomic limit.
#
# Stage C.2 will add the NLIE residual evaluator (forward operator),
# the Trotter-limit driving function `log lambda_0(x)`, and the Picard
# iteration solver. Stage C.3 will wire the production fetch dispatch.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Auxiliary-function container
# ═══════════════════════════════════════════════════════════════════════════════

"""
    JKSAuxFunctions(N)

Container for the four JKS auxiliary functions evaluated on an N-point
discretised contour:

- `b::Vector{ComplexF64}`     — b(x_j + i gamma)     (eq 27, 51)
- `b_bar::Vector{ComplexF64}` — bar b(x_j - i gamma)
- `c::Vector{ComplexF64}`     — c(x_j + i 0+)
- `c_bar::Vector{ComplexF64}` — bar c(x_j - i 0+)

Each vector has length `N`. The constructor `JKSAuxFunctions(N)` returns
zero-initialised arrays — see `init_atomic_limit!` for a physically
meaningful starting point.
"""
struct JKSAuxFunctions
    b::Vector{ComplexF64}
    b_bar::Vector{ComplexF64}
    c::Vector{ComplexF64}
    c_bar::Vector{ComplexF64}
    function JKSAuxFunctions(N::Int)
        N > 0 || throw(DomainError(N, "N must be > 0"))
        return new(
            zeros(ComplexF64, N),
            zeros(ComplexF64, N),
            zeros(ComplexF64, N),
            zeros(ComplexF64, N),
        )
    end
end

"""
    Base.length(aux::JKSAuxFunctions) -> Int

Number of contour points (equal across all four arrays).
"""
Base.length(aux::JKSAuxFunctions) = length(aux.b)

"""
    Base.copy(aux::JKSAuxFunctions) -> JKSAuxFunctions

Deep-copy of all four auxiliary-function arrays. Used by the Stage C.2
Picard iteration to keep prior-iterate values for mixing.
"""
function Base.copy(aux::JKSAuxFunctions)
    new_aux = JKSAuxFunctions(length(aux))
    copyto!(new_aux.b, aux.b)
    copyto!(new_aux.b_bar, aux.b_bar)
    copyto!(new_aux.c, aux.c)
    copyto!(new_aux.c_bar, aux.c_bar)
    return new_aux
end

# ═══════════════════════════════════════════════════════════════════════════════
# Atomic-limit initialization
# ═══════════════════════════════════════════════════════════════════════════════
#
# At t = 0 the Hubbard sites decouple. The four Boltzmann weights
#
#   z_empty  = 1
#   z_up     = exp(beta * (mu + h/2))
#   z_down   = exp(beta * (mu - h/2))
#   z_double = exp(beta * (2 mu - U))
#
# give a per-site partition function Z_site. In this limit the auxiliary
# functions become x-independent constants determined by the JKS Sec 7.3
# high-T analysis. Without the full eq (138) at hand, we use the spin-
# channel and charge-channel ratios as physically motivated initial
# guesses:
#
#   b_const    = (z_up + z_down) / (z_empty + z_double)
#   c_const    = z_up / Z_site
#   cbar_const = z_down / Z_site
#
# These are NOT claimed to be exact NLIE solutions for t > 0; they are
# the starting point Stage C.2's Picard iteration relaxes to the true
# fixed-point auxiliary functions.

"""
    init_atomic_limit!(aux, beta, U, mu; h=0.0) -> JKSAuxFunctions

Fill `aux` with x-independent constants derived from the atomic-limit
Boltzmann weights at the requested `(beta, U, mu, h)`. Returns `aux`
for chaining.

Used as the Picard iteration starting point in Stage C.2. The constants
are physically motivated but NOT verified against JKS eq (138); Stage
C.2 may revise them once the paper's high-T analysis is fully ported.
"""
function init_atomic_limit!(
    aux::JKSAuxFunctions, beta::Real, U::Real, mu::Real; h::Real=0.0
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    z_empty = 1.0
    z_up = exp(beta * (mu + h/2))
    z_down = exp(beta * (mu - h/2))
    z_double = exp(beta * (2 * mu - U))
    Z_site = z_empty + z_up + z_down + z_double

    # Paper-precise atomic-limit aux init.
    # Constraints at beta -> 0: c*cbar = exp(-betaU) and c/cbar = exp(betah),
    # so c = exp(-betaU/2 + betah/2) and cbar = exp(-betaU/2 - betah/2).
    # The b ratio (z_up+z_down)/(z_empty+z_double) is the paper-consistent leading value.
    b_const = ComplexF64((z_up + z_down) / (z_empty + z_double))
    c_const = ComplexF64(exp(-beta * U / 2 + beta * h / 2))
    cbar_const = ComplexF64(exp(-beta * U / 2 - beta * h / 2))

    fill!(aux.b, b_const)
    fill!(aux.b_bar, b_const)
    fill!(aux.c, c_const)
    fill!(aux.c_bar, cbar_const)
    return aux
end

"""
    init_atomic_limit(grid, beta, U, mu; h=0.0) -> JKSAuxFunctions

Convenience: allocate a fresh `JKSAuxFunctions(grid.N)` and initialise
it via `init_atomic_limit!`.
"""
function init_atomic_limit(grid::JKSContourGrid, beta::Real, U::Real, mu::Real; h::Real=0.0)
    aux = JKSAuxFunctions(grid.N)
    return init_atomic_limit!(aux, beta, U, mu; h=h)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.2 — elementary function phi + free-energy evaluator
#
# After the Read-tool extraction of JKS Sec 4-5 image pages, the precise
# forms are now in hand:
#
#   (1) Elementary function phi (eq 48):
#         ln phi(s) = -2 beta |s| sqrt(1 - 1/s^2)
#       for real s with |s| > 1; analytically continued elsewhere.
#       Used inside the driving terms psi_c, psi_cbar.
#
#   (2) Free energy formula (eq 49, simplest form):
#         ln Lambda = -2 pi i (U/4)
#                   + integral_L [ln z(s)]'      ln(1 + c + cbar) ds
#                   + integral_L [ln z(s - 2i eta)]' ln C(s) ds
#       with z(s) = i s (1 + sqrt(1 - 1/s^2)) (eq 23) and eta = U/4.
#       The derivative simplifies to
#         [ln z(s)]' = 1 / sqrt(s^2 - 1).
#       Contour L is just above and below the branch cut on [-1, 1];
#       the closed-contour integral picks up the discontinuity.
#
#   (3) Per-site free energy (with N -> infty Trotter limit):
#         f(T, U, mu, h) = -T ln Lambda
#
# Stage C.3 will add the NLIE residual evaluator + Picard solver +
# production fetch dispatch. The infrastructure shipped here is the
# input that solver will need at every iteration.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Elementary function phi (eq 48)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks_log_phi(s::Real, beta::Real) -> Float64

JKS elementary function on the real axis (eq 48):

    ln phi(s) = -2 beta |s| sqrt(1 - 1/s^2),    valid for |s| > 1.

For |s| <= 1 the argument of the square root is negative and the
function takes an imaginary value (branch cut). Stage C.2 returns
`-Inf` there as a sentinel — the contour integrals only touch the
branch cut along [-1, 1] and use the discontinuity rather than the
value, so the cut is never evaluated naively.
"""
function jks_log_phi(s::Real, beta::Real)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    if abs(s) <= 1.0
        return -Inf
    end
    return -2 * beta * abs(s) * sqrt(1 - 1/s^2)
end

"""
    jks_phi(s::Real, beta::Real) -> Float64

Exponential form: `phi(s) = exp(ln phi(s))`. Returns `0.0` for
`|s| <= 1` (the cut). For numerical use prefer `jks_log_phi` to avoid
underflow at large `beta`.
"""
function jks_phi(s::Real, beta::Real)
    return exp(jks_log_phi(s, beta))
end

# ═══════════════════════════════════════════════════════════════════════════════
# log z(s) derivative (eq 23 simplified)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks_log_z_deriv(s::Number) -> ComplexF64

Derivative `[ln z(s)]' = 1 / sqrt(s^2 - 1)`, where `z(s) = i s (1 + sqrt(1 - 1/s^2))`
(eq 23). On the real axis `|s| > 1` this is real; for `|s| < 1` it is
imaginary (branch cut along [-1, 1]).

This is the kernel of the closed-contour integration in the JKS free
energy formula eq (49).
"""
function jks_log_z_deriv(s::Number)
    return 1 / sqrt(s^2 - 1 + 0im)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Free energy evaluator (eq 49)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    free_energy_jks(aux, grid, beta, U; mu=U/2) -> Float64

Per-site free energy of the 1D Hubbard model from JKS eq (49), given
the auxiliary functions `aux` on the discretized contour `grid`.

This evaluator is correct **once aux is the converged NLIE solution**;
on un-converged aux (such as Stage C.1's atomic-limit initial guess)
it returns the free energy "as if aux were the answer" — a useful
baseline for atomic-limit consistency checks, but not the physical f
for t > 0.

# Arguments

- `aux`: JKSAuxFunctions container holding (b, b_bar, c, c_bar) sampled
  on `grid.x`.
- `grid`: JKSContourGrid. **Note**: in Stage C.2 we interpret
  `grid.gamma` as the Hubbard parameter `eta = U / 4` (the kernel-pole
  shift), not as a free contour parameter. The independent contour
  shift parameter `alpha` is not used here (it only enters the NLIE
  solver of Stage C.3).
- `beta`, `U`, `mu`: physics parameters. Default `mu = U/2` is the
  half-filling case.

# Returns

The per-site free-energy density `f = -T ln Lambda / (2 pi i)`,
evaluated by the discrete Simpson-style quadrature of eq (49) over
the branch cut `[-1, 1]`.

# References

- JKS 1998 eq (23), (38), (49) — z(s), kernel, free energy
"""
function free_energy_jks(
    aux::JKSAuxFunctions, grid::JKSContourGrid, beta::Real, U::Real; mu::Real=U/2
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    length(aux) == grid.N ||
        throw(DimensionMismatch("aux length $(length(aux)) != grid.N $(grid.N)"))

    eta = U / 4

    # Restrict to grid points inside the branch cut [-1, 1]:
    # the contour integral only contributes on the cut.
    cut_mask = [-1.0 <= x <= 1.0 for x in grid.x]

    # First integral:  -2i int_{-1}^{1} ln(1 + c + cbar) / sqrt(1 - s^2) ds
    integrand_1 = ComplexF64[
        if cut_mask[j]
            log((1 + aux.c[j] + aux.c_bar[j]) / aux.c[j]) / sqrt(1 - grid.x[j]^2)
        else
            0.0 + 0im
        end for j in 1:grid.N
    ]
    int_1 = -2im * sum(integrand_1) * grid.dx

    # Second integral:  oint [ln z(s - 2 i eta)]' ln C(s) ds.
    # The pole at s - 2i eta = ±1 is off the real axis for eta > 0; this
    # second integral is regular on [-1, 1] and we evaluate it as a
    # straight Riemann sum.
    integrand_2 = ComplexF64[
        if cut_mask[j]
            jks_log_z_deriv(grid.x[j] - 2im * eta) * log(1 + aux.c[j])
        else
            0.0 + 0im
        end for j in 1:grid.N
    ]
    int_2 = sum(integrand_2) * grid.dx

    # ln Lambda from eq (49):
    ln_Lambda = -2pi * im * (U/4) + int_1 + int_2

    # Per-site f = -T ln Lambda / (2 pi i)  — the factor 2 pi i in
    # the LHS of eq (49) means ln Lambda on the RHS is unnormalised;
    # divide by 2 pi i to get the physical eigenvalue.
    f = (1/beta) * ln_Lambda / (2pi * im)  # sign per JKS convention (Stage C.10 fix)

    return real(f)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.3 partial — complex analytic continuation of phi
#
# JKS Sec 5 uses the elementary function phi evaluated on contours shifted
# off the real axis by +- i alpha (alpha is the b/b_bar contour shift, free
# with 0 < alpha < eta = U/4). The Stage C.2 jks_log_phi only covered the
# real axis at |s| > 1; this file adds the complex analytic continuation
#
#     ln phi(s) = -2 beta s sqrt(1 - 1/s^2)
#
# with the principal branch of sqrt (positive real part for large |s|).
#
# This is the input for the driving terms psi_b, psi_c, psi_c_bar (eq 58-59)
# which Stage C.4 will assemble.
# ─────────────────────────────────────────────────────────────────────────────

"""
    jks_log_phi_complex(s::Number, beta::Real) -> ComplexF64

Complex analytic continuation of the JKS elementary function phi (eq 48):

    ln phi(s) = -2 beta s sqrt(1 - 1/s^2)

For real `s` with |s| > 1 this reduces to the real-axis form
`ln phi(s) = -2 beta |s| sqrt(1 - 1/s^2)`. For real `s` with |s| <= 1
the value is pure imaginary — this is the branch-cut discontinuity
used by the contour integrals in eq (49).

For complex `s = x + i y` with `y != 0`, returns the value on the
principal sheet.

Throws `DomainError` for `beta <= 0`.
"""
function jks_log_phi_complex(s::Number, beta::Real)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    return -2 * beta * sqrt(s^2 - 1 + 0im)
end

"""
    jks_phi_complex(s::Number, beta::Real) -> ComplexF64

Exponential form: `phi(s) = exp(ln phi(s))` with complex `s`. The result
is finite for any `s` (the cut is a discontinuity, not a singularity).
"""
function jks_phi_complex(s::Number, beta::Real)
    return exp(jks_log_phi_complex(s, beta))
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.4 — driving terms + NLIE residual + Picard solver
#
# Closes the bulk implementation of the JKS Hubbard finite-T NLIE except
# for the production dispatch switch and Fig 5-10 agreement tests
# (deferred to Stage C.5).
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Driving terms (eq 58-59)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks_driving_b(grid, beta, U, alpha; H=0.0) -> Vector{ComplexF64}

JKS driving function `psi_b(x_j) = -beta U - beta H + log phi(x+i alpha) - log phi(x-i alpha)` (eq 59).
"""
function jks_driving_b(grid::JKSContourGrid, beta::Real, U::Real, alpha::Real; H::Real=0.0)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    alpha > 0 || throw(DomainError(alpha, "alpha must be > 0"))
    return [
        -beta * U - beta * H + jks_log_phi_complex(x + im * alpha, beta) -
        jks_log_phi_complex(x - im * alpha, beta) for x in grid.x
    ]
end

"""
    jks_driving_c(grid, beta, U, mu; H=0.0) -> Vector{ComplexF64}

JKS driving function `psi_c(x_j) = -beta U/2 - beta(mu + H/2) + log phi(x)` (eq 58).
"""
function jks_driving_c(grid::JKSContourGrid, beta::Real, U::Real, mu::Real; H::Real=0.0)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    return [
        -beta * U / 2 + beta * (mu + H/2) + jks_log_phi_complex(x + 0im, beta) for
        x in grid.x
    ]
end

"""
    jks_driving_cbar(grid, beta, U, mu; H=0.0) -> Vector{ComplexF64}

Particle-hole conjugate of `jks_driving_c`. Satisfies `psi_c + psi_cbar = -beta U`.
"""
function jks_driving_cbar(grid::JKSContourGrid, beta::Real, U::Real, mu::Real; H::Real=0.0)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    return [
        -beta * U / 2 - beta * (mu + H/2) - jks_log_phi_complex(x + 0im, beta) for
        x in grid.x
    ]
end

# ═══════════════════════════════════════════════════════════════════════════════
# NLIE residual evaluator (b-channel of eq 53)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks_nlie_residual(aux, grid, beta, U, mu, alpha; H=0.0) -> Vector{ComplexF64}

Compute the b-channel NLIE residual `log b - RHS` for the given `aux`.
The L-infinity norm of the residual measures convergence; at the
converged NLIE solution it is zero up to discretization error.

This is the b-equation residual only — Stage C.5 will extend to the
full (b, c, c_bar) triple.
"""
function jks_nlie_residual(
    aux::JKSAuxFunctions,
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real,
    alpha::Real;
    H::Real=0.0,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    alpha > 0 || throw(DomainError(alpha, "alpha must be > 0"))
    length(aux) == grid.N ||
        throw(DimensionMismatch("aux length $(length(aux)) != grid.N $(grid.N)"))

    K1 = build_kernel_matrix(grid, 1)
    K2 = build_kernel_matrix(grid, 2)

    psi_b = jks_driving_b(grid, beta, U, alpha; H=H)

    log_B = log.(1 .+ aux.b)
    log_Bbar = log.(1 .+ 1 ./ aux.b_bar)
    log_c = log.(aux.c)
    log_C = log.(1 .+ aux.c)
    log_c_minus_C = log_c .- log_C

    # Paper eq (47) b residual: log b = -betaH + K_2 ⊛ log B + K_1 ⊛ (log c - log C)
    rhs = (-beta * H) .+ apply_kernel(K2, log_B) + apply_kernel(K1, log_c_minus_C)

    log_b = log.(aux.b)
    return log_b .- rhs
end

# ═══════════════════════════════════════════════════════════════════════════════
# Picard solver (b-channel only, alpha-mixing)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    JKSSolution

Container for `solve_jks_nlie_b_only` result.
"""
struct JKSSolution
    aux::JKSAuxFunctions
    iterations::Int
    residual::Float64
    converged::Bool
end

"""
    solve_jks_nlie_b_only(grid, beta, U, mu; alpha=U/6, H=0, tol=1e-6,
                          maxiter=200, alpha_mix=0.3) -> JKSSolution

Picard iteration on the b-channel of the JKS NLIE only (Stage C.4
scope). c and c_bar are held at the atomic-limit constants. Init from
`init_atomic_limit`; update rule

    log b_new = log b_old - alpha_mix * residual

until L-infinity norm of residual is below `tol`.
"""
function solve_jks_nlie_b_only(
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real;
    alpha::Real=U/6,
    H::Real=0.0,
    tol::Real=1e-6,
    maxiter::Int=200,
    alpha_mix::Real=0.3,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    0 < alpha_mix <= 1 || throw(DomainError(alpha_mix, "alpha_mix must be in (0, 1]"))
    aux = init_atomic_limit(grid, beta, U, mu; h=H)

    last_residual = Inf
    for iter in 1:maxiter
        res = jks_nlie_residual(aux, grid, beta, U, mu, alpha; H=H)
        last_residual = maximum(abs.(res))
        if last_residual < tol
            return JKSSolution(aux, iter, last_residual, true)
        end
        log_b_new = log.(aux.b) .- alpha_mix .* res
        aux.b .= exp.(log_b_new)
        aux.b_bar .= aux.b
    end
    return JKSSolution(aux, maxiter, last_residual, false)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.5 — physical kernel regularization + improved solver
#
# Stage C.4 build_kernel_matrix used eps = 1e-10 as the diagonal pole
# regularizer, giving K[j,j] ~ 4e8 which dominated the residual. The
# physically correct regularization is the contour shift alpha (with
# 0 < alpha < eta), giving K(0 + i alpha) of order O(1).
# ─────────────────────────────────────────────────────────────────────────────

"""
    build_kernel_matrix_shifted(grid, n, alpha_shift) -> Matrix{ComplexF64}

Kernel matrix with diagonal regularized by the physical contour shift
`alpha_shift` (typically 0 < alpha_shift < eta = U/4) rather than the
tiny `eps = 1e-10` used by `build_kernel_matrix`.

K[j, k] = jks_kernel_K_n_concrete(x_j - x_k + i * alpha_shift, n, gamma) * dx

The j = k diagonal evaluates at K(i * alpha_shift), which is O(1) for
alpha_shift ~ eta.
"""
function build_kernel_matrix_shifted(grid::JKSContourGrid, n::Integer, alpha_shift::Real)
    n > 0 || throw(DomainError(n, "kernel index n must be > 0"))
    alpha_shift > 0 || throw(DomainError(alpha_shift, "alpha_shift must be > 0"))
    N = grid.N
    K = zeros(ComplexF64, N, N)
    for j in 1:N, k in 1:N
        K[j, k] =
            jks_kernel_K_n_concrete(
                grid.x[j] - grid.x[k] + im * alpha_shift, n, grid.gamma
            ) * grid.dx
    end
    return K
end

# Stage C.18: K_bar kernel with -i*alpha shift (conjugate of K_n).
function build_kernel_matrix_shifted_bar(
    grid::JKSContourGrid, n::Integer, alpha_shift::Real
)
    n > 0 || throw(DomainError(n, "kernel index n must be > 0"))
    alpha_shift > 0 || throw(DomainError(alpha_shift, "alpha_shift must be > 0"))
    N = grid.N
    K = zeros(ComplexF64, N, N)
    for j in 1:N, k in 1:N
        K[j, k] =
            jks_kernel_K_n_concrete(
                grid.x[j] - grid.x[k] - im * alpha_shift, n, grid.gamma
            ) * grid.dx
    end
    return K
end

"""
    jks_nlie_residual_shifted(aux, grid, beta, U, mu, alpha; H=0.0) -> Vector{ComplexF64}

Same as `jks_nlie_residual` but uses the physically-regularized kernels
via `build_kernel_matrix_shifted` with the same `alpha` for both the
contour shift and the kernel regularizer.
"""
function jks_nlie_residual_shifted(
    aux::JKSAuxFunctions,
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real,
    alpha::Real;
    H::Real=0.0,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    alpha > 0 || throw(DomainError(alpha, "alpha must be > 0"))
    length(aux) == grid.N ||
        throw(DimensionMismatch("aux length $(length(aux)) != grid.N $(grid.N)"))

    K1 = build_kernel_matrix_shifted(grid, 1, alpha)
    K2 = build_kernel_matrix_shifted(grid, 2, alpha)

    psi_b = jks_driving_b(grid, beta, U, alpha; H=H)

    log_B = log.(1 .+ aux.b)
    log_Bbar = log.(1 .+ 1 ./ aux.b_bar)
    log_c = log.(aux.c)
    log_C = log.(1 .+ aux.c)
    log_c_minus_C = log_c .- log_C

    # Paper eq (47) b residual: log b = -betaH + K_2 ⊛ log B + K_1 ⊛ (log c - log C)
    rhs = (-beta * H) .+ apply_kernel(K2, log_B) + apply_kernel(K1, log_c_minus_C)

    log_b = log.(aux.b)
    return log_b .- rhs
end

"""
    solve_jks_nlie_shifted(grid, beta, U, mu; alpha=U/6, H=0, tol=1e-6,
                          maxiter=200, alpha_mix=0.05) -> JKSSolution

Picard solver using the physically-shifted kernels. Default `alpha_mix`
is 0.05 (smaller than Stage C.4's 0.3) because even with correct kernel
regularization the Jacobian can drive Picard unstable for large mixing.
"""
function solve_jks_nlie_shifted(
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real;
    alpha::Real=U/6,
    H::Real=0.0,
    tol::Real=1e-6,
    maxiter::Int=200,
    alpha_mix::Real=0.05,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    0 < alpha_mix <= 1 || throw(DomainError(alpha_mix, "alpha_mix must be in (0, 1]"))
    aux = init_atomic_limit(grid, beta, U, mu; h=H)

    last_residual = Inf
    for iter in 1:maxiter
        res = jks_nlie_residual_shifted(aux, grid, beta, U, mu, alpha; H=H)
        last_residual = maximum(abs.(res))
        if !isfinite(last_residual)
            return JKSSolution(aux, iter, last_residual, false)
        end
        if last_residual < tol
            return JKSSolution(aux, iter, last_residual, true)
        end
        log_b_new = log.(aux.b) .- alpha_mix .* res
        aux.b .= exp.(log_b_new)
        aux.b_bar .= aux.b
    end
    return JKSSolution(aux, maxiter, last_residual, false)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.6 — adaptive alpha-mix Picard solver
#
# Stage C.5 fixed the kernel scale; adaptive mixing extends convergence
# from the high-T regime down into intermediate T / U.
# ─────────────────────────────────────────────────────────────────────────────

"""
    solve_jks_nlie_adaptive(grid, beta, U, mu; alpha=U/6, H=0, tol=1e-6,
                            maxiter=2000, alpha_mix_init=0.02,
                            alpha_mix_floor=1e-6, alpha_mix_cap=0.5,
                            shrink=0.5, grow=1.1) -> JKSSolution

Adaptive alpha-mix Picard solver. Updates the mixing parameter based on
whether the residual is decreasing:

  - shrink alpha_mix by `shrink` if residual grows by > 10%
  - grow alpha_mix by `grow` if residual shrinks by > 10%
  - floor at `alpha_mix_floor`
  - cap at `alpha_mix_cap`
"""
function solve_jks_nlie_adaptive(
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real;
    alpha::Real=U/6,
    H::Real=0.0,
    tol::Real=1e-6,
    maxiter::Int=2000,
    alpha_mix_init::Real=0.02,
    alpha_mix_floor::Real=1e-6,
    alpha_mix_cap::Real=0.5,
    shrink::Real=0.5,
    grow::Real=1.1,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    0 < alpha_mix_init <= 1 ||
        throw(DomainError(alpha_mix_init, "alpha_mix_init must be in (0, 1]"))
    0 < shrink < 1 || throw(DomainError(shrink, "shrink must be in (0, 1)"))
    grow > 1 || throw(DomainError(grow, "grow must be > 1"))

    aux = init_atomic_limit(grid, beta, U, mu; h=H)

    alpha_mix = alpha_mix_init
    prev_residual = Inf
    last_residual = Inf

    for iter in 1:maxiter
        res = jks_nlie_residual_shifted(aux, grid, beta, U, mu, alpha; H=H)
        last_residual = maximum(abs.(res))

        if !isfinite(last_residual)
            return JKSSolution(aux, iter, last_residual, false)
        end
        if last_residual < tol
            return JKSSolution(aux, iter, last_residual, true)
        end

        if last_residual > prev_residual * 1.1
            alpha_mix = max(alpha_mix * shrink, alpha_mix_floor)
        elseif last_residual < prev_residual * 0.9
            alpha_mix = min(alpha_mix * grow, alpha_mix_cap)
        end

        log_b_new = log.(aux.b) .- alpha_mix .* res
        aux.b .= exp.(log_b_new)
        aux.b_bar .= aux.b

        prev_residual = last_residual
    end
    return JKSSolution(aux, maxiter, last_residual, false)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.7 — damped Newton solver
#
# Picard variants (Stage C.4, C.5, C.6) are first-order. Newton gives
# quadratic convergence near the solution at the price of a Jacobian
# solve per iteration. We use forward finite differences for the
# Jacobian and damped line search to handle overshoot.
# ─────────────────────────────────────────────────────────────────────────────

"""
    jks_jacobian_b_finite_diff(aux, grid, beta, U, mu, alpha; H=0.0, eps=1e-6)
        -> Matrix{ComplexF64}

Numerical Jacobian `J[i, k] = d residual_i / d log b_k` via forward
finite difference at step `eps`. Requires `N + 1` residual evaluations
on an `N`-point grid.
"""
function jks_jacobian_b_finite_diff(
    aux::JKSAuxFunctions,
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real,
    alpha::Real;
    H::Real=0.0,
    eps::Real=1e-6,
)
    N = grid.N
    res0 = jks_nlie_residual_shifted(aux, grid, beta, U, mu, alpha; H=H)
    J = zeros(ComplexF64, N, N)
    log_b0 = log.(aux.b)
    for k in 1:N
        aux_pert = copy(aux)
        log_b_pert = copy(log_b0)
        log_b_pert[k] += eps
        aux_pert.b .= exp.(log_b_pert)
        aux_pert.b_bar .= aux_pert.b
        res_pert = jks_nlie_residual_shifted(aux_pert, grid, beta, U, mu, alpha; H=H)
        J[:, k] = (res_pert .- res0) ./ eps
    end
    return J
end

"""
    solve_jks_nlie_newton(grid, beta, U, mu; alpha=U/6, H=0, tol=1e-6,
                          maxiter=50, jac_eps=1e-6,
                          damps=(1.0, 0.5, 0.25, 0.1, 0.05, 0.01)) -> JKSSolution

Damped Newton solver. At each iteration: compute residual r and Jacobian
J, solve `J · δ = -r`, then try damping factors from largest to smallest
and accept the first that decreases the residual.
"""
function solve_jks_nlie_newton(
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real;
    alpha::Real=U/6,
    H::Real=0.0,
    tol::Real=1e-6,
    maxiter::Int=50,
    jac_eps::Real=1e-6,
    damps=(1.0, 0.5, 0.25, 0.1, 0.05, 0.01),
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    aux = init_atomic_limit(grid, beta, U, mu; h=H)
    last_residual = Inf

    for iter in 1:maxiter
        res = jks_nlie_residual_shifted(aux, grid, beta, U, mu, alpha; H=H)
        last_residual = maximum(abs.(res))

        if !isfinite(last_residual)
            return JKSSolution(aux, iter, last_residual, false)
        end
        if last_residual < tol
            return JKSSolution(aux, iter, last_residual, true)
        end

        J = jks_jacobian_b_finite_diff(aux, grid, beta, U, mu, alpha; H=H, eps=jac_eps)

        delta = try
            J \ (-res)
        catch
            return JKSSolution(aux, iter, last_residual, false)
        end

        if !all(isfinite, delta)
            return JKSSolution(aux, iter, last_residual, false)
        end

        log_b_old = log.(aux.b)
        accepted = false
        for damp in damps
            log_b_new = log_b_old .+ damp .* delta
            aux_try = copy(aux)
            aux_try.b .= exp.(log_b_new)
            aux_try.b_bar .= aux_try.b
            res_try = try
                jks_nlie_residual_shifted(aux_try, grid, beta, U, mu, alpha; H=H)
            catch
                continue
            end
            res_try_norm = maximum(abs.(res_try))
            if isfinite(res_try_norm) && res_try_norm < last_residual
                aux.b .= aux_try.b
                aux.b_bar .= aux_try.b_bar
                accepted = true
                break
            end
        end
        if !accepted
            return JKSSolution(aux, iter, last_residual, false)
        end
    end
    return JKSSolution(aux, maxiter, last_residual, false)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.8 — beta-continuation Newton solver
#
# Stage C.7 showed Newton converges in 2 iter at beta = 0.01 but stalls
# at beta = 1.0 because the atomic-limit init becomes a bad guess.
# Continuation: solve at high T, then walk beta up using each converged
# aux as the next initial guess.
# ─────────────────────────────────────────────────────────────────────────────

"""
    solve_jks_nlie_newton_from(aux_init, grid, beta, U, mu; alpha, H, tol,
                                maxiter, jac_eps, damps) -> JKSSolution

Like `solve_jks_nlie_newton` but takes an explicit initial guess
`aux_init` instead of allocating from `init_atomic_limit`. Used as the
inner loop of the beta-continuation solver.
"""
function solve_jks_nlie_newton_from(
    aux_init::JKSAuxFunctions,
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real;
    alpha::Real=U/6,
    H::Real=0.0,
    tol::Real=1e-6,
    maxiter::Int=50,
    jac_eps::Real=1e-6,
    damps=(1.0, 0.5, 0.25, 0.1, 0.05, 0.01),
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    aux = copy(aux_init)
    last_residual = Inf

    for iter in 1:maxiter
        res = jks_nlie_residual_shifted(aux, grid, beta, U, mu, alpha; H=H)
        last_residual = maximum(abs.(res))
        if !isfinite(last_residual)
            return JKSSolution(aux, iter, last_residual, false)
        end
        if last_residual < tol
            return JKSSolution(aux, iter, last_residual, true)
        end

        J = jks_jacobian_b_finite_diff(aux, grid, beta, U, mu, alpha; H=H, eps=jac_eps)
        delta = try
            J \ (-res)
        catch
            return JKSSolution(aux, iter, last_residual, false)
        end
        if !all(isfinite, delta)
            return JKSSolution(aux, iter, last_residual, false)
        end

        log_b_old = log.(aux.b)
        accepted = false
        for damp in damps
            log_b_new = log_b_old .+ damp .* delta
            aux_try = copy(aux)
            aux_try.b .= exp.(log_b_new)
            aux_try.b_bar .= aux_try.b
            res_try = try
                jks_nlie_residual_shifted(aux_try, grid, beta, U, mu, alpha; H=H)
            catch
                continue
            end
            res_try_norm = maximum(abs.(res_try))
            if isfinite(res_try_norm) && res_try_norm < last_residual
                aux.b .= aux_try.b
                aux.b_bar .= aux_try.b_bar
                accepted = true
                break
            end
        end
        if !accepted
            return JKSSolution(aux, iter, last_residual, false)
        end
    end
    return JKSSolution(aux, maxiter, last_residual, false)
end

"""
    solve_jks_nlie_continuation(grid, beta_target, U, mu; alpha=U/6, H=0,
                                 beta_start=0.01, step_init=0.5,
                                 step_min=1e-3, step_max=1.0,
                                 grow_factor=1.2, shrink_factor=0.5,
                                 tol=1e-6, inner_maxiter=30,
                                 outer_maxsteps=200) -> JKSSolution

Beta-continuation Newton solver. Solve at `beta_start` from atomic-limit
init, then walk `beta` up multiplicatively (`beta_new = beta * (1 + step)`)
to `beta_target`. On inner-Newton convergence: grow `step`. On failure:
shrink `step`. Give up if `step < step_min`.
"""
function solve_jks_nlie_continuation(
    grid::JKSContourGrid,
    beta_target::Real,
    U::Real,
    mu::Real;
    alpha::Real=U/6,
    H::Real=0.0,
    beta_start::Real=0.01,
    step_init::Real=0.5,
    step_min::Real=1e-3,
    step_max::Real=1.0,
    grow_factor::Real=1.2,
    shrink_factor::Real=0.5,
    tol::Real=1e-6,
    inner_maxiter::Int=30,
    outer_maxsteps::Int=200,
)
    beta_target > 0 || throw(DomainError(beta_target, "beta_target must be > 0"))
    beta_start > 0 || throw(DomainError(beta_start, "beta_start must be > 0"))
    beta_start <= beta_target || throw(ArgumentError("beta_start must be <= beta_target"))
    0 < shrink_factor < 1 ||
        throw(DomainError(shrink_factor, "shrink_factor must be in (0, 1)"))
    grow_factor > 1 || throw(DomainError(grow_factor, "grow_factor must be > 1"))

    aux = init_atomic_limit(grid, beta_start, U, mu; h=H)
    sol_init = solve_jks_nlie_newton_from(
        aux, grid, beta_start, U, mu; alpha=alpha, H=H, tol=tol, maxiter=inner_maxiter
    )
    if !sol_init.converged
        return JKSSolution(sol_init.aux, sol_init.iterations, sol_init.residual, false)
    end
    aux = sol_init.aux

    beta_current = beta_start
    step = step_init
    total_iter = sol_init.iterations

    for outer_iter in 1:outer_maxsteps
        if beta_current >= beta_target - 1e-12
            return JKSSolution(aux, total_iter, sol_init.residual, true)
        end

        beta_try = min(beta_current * (1 + step), beta_target)
        sol = solve_jks_nlie_newton_from(
            aux, grid, beta_try, U, mu; alpha=alpha, H=H, tol=tol, maxiter=inner_maxiter
        )
        total_iter += sol.iterations

        if sol.converged
            beta_current = beta_try
            aux = sol.aux
            step = min(step * grow_factor, step_max)
        else
            step *= shrink_factor
            if step < step_min
                return JKSSolution(aux, total_iter, sol.residual, false)
            end
        end
    end

    final_res = maximum(
        abs.(jks_nlie_residual_shifted(aux, grid, beta_current, U, mu, alpha; H=H))
    )
    return JKSSolution(aux, total_iter, final_res, beta_current >= beta_target - 1e-12)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.9 — high-level free energy wrapper
#
# Wraps grid construction + continuation solver + free-energy evaluator
# into a single call. Production fetch dispatch wiring (modifying
# Hubbard1D_thermal_stopgap) is deferred to Stage C.10.
#
# Limitation: JKS NLIE is derived with t = 1 normalization. Callers
# wanting general t can rescale beta and U externally.
# ─────────────────────────────────────────────────────────────────────────────

"""
    hubbard1d_jks_free_energy(t, U, mu, beta; alpha=U/6, H=0.0,
                              grid_N=16, x_max=2.0,
                              tol=1e-3, beta_start=0.01,
                              inner_maxiter=30, outer_maxsteps=200)
        -> Float64

Per-site Helmholtz free-energy density of the 1D Hubbard chain via the
JKS QTM NLIE + beta-continuation Newton solver.

Returns the free-energy density on convergence; returns `NaN` if the
continuation cannot reach `beta` (the solver stalls).

# Arguments

- `t = 1`: enforced. The JKS NLIE assumes the t = 1 normalization;
  callers wanting general t pass U/t and mu/t, multiply beta by t.
- `U`, `mu`: Hubbard parameters; half-filling is mu = U/2.
- `beta`: inverse temperature.
- `alpha`: contour shift, free in (0, U/4).
- `H`: magnetic field.
- `grid_N`, `x_max`: discretization grid size + half-width.
- `tol`, `beta_start`, `inner_maxiter`, `outer_maxsteps`: continuation knobs.
"""
function hubbard1d_jks_free_energy(
    t::Real,
    U::Real,
    mu::Real,
    beta::Real;
    alpha::Real=U/6,
    H::Real=0.0,
    grid_N::Int=16,
    x_max::Real=2.0,
    tol::Real=1e-3,
    beta_start::Real=0.01,
    inner_maxiter::Int=30,
    outer_maxsteps::Int=200,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    U >= 0 || throw(DomainError(U, "U must be >= 0"))
    isapprox(t, 1.0) ||
        throw(ArgumentError("JKS wrapper assumes t = 1; rescale externally for general t"))

    eta = U / 4

    if alpha >= eta
        return NaN
    end

    grid = JKSContourGrid(grid_N, eta; x_max=x_max)

    bs = min(beta_start, beta)
    sol = solve_jks_nlie_continuation(
        grid,
        beta,
        U,
        mu;
        alpha=alpha,
        H=H,
        beta_start=bs,
        tol=tol,
        inner_maxiter=inner_maxiter,
        outer_maxsteps=outer_maxsteps,
    )

    if !sol.converged
        return NaN
    end

    return free_energy_jks(sol.aux, grid, beta, U; mu=mu)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.11 — c, c_bar channel residuals + 3-equation system
#
# Stage C.4-C.10 only iterated the b channel; c and c_bar were held at
# atomic-limit constants. That left a ~3 % offset in the high-T
# hubbard1d_jks_free_energy vs atomic_free_energy comparison (#544).
#
# Adds residuals for the c and c_bar equations from JKS eq (53) and a
# concatenation `jks_nlie_residual_full` that Stage C.12 will drive to
# zero via a 3N x 3N Newton solver.
# ─────────────────────────────────────────────────────────────────────────────

"""
    jks_nlie_residual_c(aux, grid, beta, U, mu, alpha; H=0.0) -> Vector{ComplexF64}

c-channel NLIE residual `log c+ - RHS_c` per JKS eq (53).
"""
function jks_nlie_residual_c(
    aux::JKSAuxFunctions,
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real,
    alpha::Real;
    H::Real=0.0,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    alpha > 0 || throw(DomainError(alpha, "alpha must be > 0"))
    length(aux) == grid.N ||
        throw(DimensionMismatch("aux length $(length(aux)) != grid.N $(grid.N)"))

    # Paper eq (47): log c = psi_c - K_1 ⊓⊔ log B - K_1 ◦ log C
    K1 = build_kernel_matrix_shifted(grid, 1, alpha)

    psi_c = jks_driving_c(grid, beta, U, mu; H=H)

    log_B = log.(1 .+ aux.b)
    log_C = log.(1 .+ aux.c)

    rhs = psi_c - apply_kernel(K1, log_B) - apply_kernel(K1, log_C)

    log_c = log.(aux.c)
    return log_c .- rhs
end

"""
    jks_nlie_residual_cbar(aux, grid, beta, U, mu, alpha; H=0.0) -> Vector{ComplexF64}

c_bar-channel NLIE residual (particle-hole conjugate structure).
"""
function jks_nlie_residual_cbar(
    aux::JKSAuxFunctions,
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real,
    alpha::Real;
    H::Real=0.0,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    alpha > 0 || throw(DomainError(alpha, "alpha must be > 0"))
    length(aux) == grid.N ||
        throw(DimensionMismatch("aux length $(length(aux)) != grid.N $(grid.N)"))

    # Paper eq (47): log cbar = psi_cbar + K_1 ⊓⊔ log B + K_1 ◦ log C
    K1 = build_kernel_matrix_shifted(grid, 1, alpha)

    psi_cbar = jks_driving_cbar(grid, beta, U, mu; H=H)

    log_B = log.(1 .+ aux.b)
    log_C = log.(1 .+ aux.c)

    rhs = psi_cbar + apply_kernel(K1, log_B) + apply_kernel(K1, log_C)

    log_cbar = log.(aux.c_bar)
    return log_cbar .- rhs
end

"""
    jks_nlie_residual_full(aux, grid, beta, U, mu, alpha; H=0.0) -> Vector{ComplexF64}

Concatenated NLIE residual covering all three channels (b, c, c_bar) as
a single length-`3N` complex vector ordered `[res_b; res_c; res_cbar]`.
Stage C.12's full-system Newton solver drives this to zero.
"""
function jks_nlie_residual_full(
    aux::JKSAuxFunctions,
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real,
    alpha::Real;
    H::Real=0.0,
)
    res_b = jks_nlie_residual_shifted(aux, grid, beta, U, mu, alpha; H=H)
    res_c = jks_nlie_residual_c(aux, grid, beta, U, mu, alpha; H=H)
    res_cbar = jks_nlie_residual_cbar(aux, grid, beta, U, mu, alpha; H=H)
    return vcat(res_b, res_c, res_cbar)
end

# ─────────────────────────────────────────────────────────────────────────────
# Hubbard1DJKSNLIE Stage C.12 — full 3N x 3N Newton solver
#
# Stage C.7 b-only Newton converged in 2 iter at high T but stalled at
# mid-T because c, c_bar were held at atomic-limit constants. Stage C.11
# added the c and c_bar residuals. This file assembles them into a
# 3N x 3N Jacobian and damped Newton solver.
# ─────────────────────────────────────────────────────────────────────────────

"""
    jks_jacobian_full_finite_diff(aux, grid, beta, U, mu, alpha; H=0.0, eps=1e-6)
        -> Matrix{ComplexF64}

Numerical Jacobian of the full 3N residual (`[res_b; res_c; res_cbar]`)
with respect to `[log b; log c; log c_bar]`, via forward finite
differences. Returns a `3N x 3N` complex matrix.
"""
function jks_jacobian_full_finite_diff(
    aux::JKSAuxFunctions,
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real,
    alpha::Real;
    H::Real=0.0,
    eps::Real=1e-6,
)
    N = grid.N
    M = 3 * N
    res0 = jks_nlie_residual_full(aux, grid, beta, U, mu, alpha; H=H)
    J = zeros(ComplexF64, M, M)

    log_b0 = log.(aux.b)
    for k in 1:N
        aux_pert = copy(aux)
        log_pert = copy(log_b0)
        log_pert[k] += eps
        aux_pert.b .= exp.(log_pert)
        aux_pert.b_bar .= aux_pert.b
        res_pert = jks_nlie_residual_full(aux_pert, grid, beta, U, mu, alpha; H=H)
        J[:, k] = (res_pert .- res0) ./ eps
    end

    log_c0 = log.(aux.c)
    for k in 1:N
        aux_pert = copy(aux)
        log_pert = copy(log_c0)
        log_pert[k] += eps
        aux_pert.c .= exp.(log_pert)
        res_pert = jks_nlie_residual_full(aux_pert, grid, beta, U, mu, alpha; H=H)
        J[:, N + k] = (res_pert .- res0) ./ eps
    end

    log_cbar0 = log.(aux.c_bar)
    for k in 1:N
        aux_pert = copy(aux)
        log_pert = copy(log_cbar0)
        log_pert[k] += eps
        aux_pert.c_bar .= exp.(log_pert)
        res_pert = jks_nlie_residual_full(aux_pert, grid, beta, U, mu, alpha; H=H)
        J[:, 2 * N + k] = (res_pert .- res0) ./ eps
    end

    return J
end

"""
    solve_jks_nlie_full_newton(grid, beta, U, mu; alpha=U/6, H=0, tol=1e-6,
                                maxiter=50, jac_eps=1e-6,
                                damps=(1.0, 0.5, 0.25, 0.1, 0.05, 0.01))
        -> JKSSolution

Damped 3N x 3N Newton solver on the full (b, c, c_bar) residual system.
"""
function solve_jks_nlie_full_newton(
    grid::JKSContourGrid,
    beta::Real,
    U::Real,
    mu::Real;
    alpha::Real=U/6,
    H::Real=0.0,
    tol::Real=1e-6,
    maxiter::Int=50,
    jac_eps::Real=1e-6,
    damps=(1.0, 0.5, 0.25, 0.1, 0.05, 0.01),
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    aux = init_atomic_limit(grid, beta, U, mu; h=H)
    N = grid.N
    last_residual = Inf

    for iter in 1:maxiter
        res = jks_nlie_residual_full(aux, grid, beta, U, mu, alpha; H=H)
        last_residual = maximum(abs.(res))

        if !isfinite(last_residual)
            return JKSSolution(aux, iter, last_residual, false)
        end
        if last_residual < tol
            return JKSSolution(aux, iter, last_residual, true)
        end

        J = jks_jacobian_full_finite_diff(aux, grid, beta, U, mu, alpha; H=H, eps=jac_eps)

        delta = try
            J \ (-res)
        catch
            return JKSSolution(aux, iter, last_residual, false)
        end
        if !all(isfinite, delta)
            return JKSSolution(aux, iter, last_residual, false)
        end

        delta_b = delta[1:N]
        delta_c = delta[(N + 1):(2 * N)]
        delta_cbar = delta[(2 * N + 1):(3 * N)]

        log_b_old = log.(aux.b)
        log_c_old = log.(aux.c)
        log_cbar_old = log.(aux.c_bar)

        accepted = false
        for damp in damps
            aux_try = copy(aux)
            aux_try.b .= exp.(log_b_old .+ damp .* delta_b)
            aux_try.b_bar .= aux_try.b
            aux_try.c .= exp.(log_c_old .+ damp .* delta_c)
            aux_try.c_bar .= exp.(log_cbar_old .+ damp .* delta_cbar)
            res_try = try
                jks_nlie_residual_full(aux_try, grid, beta, U, mu, alpha; H=H)
            catch
                continue
            end
            res_try_norm = maximum(abs.(res_try))
            if isfinite(res_try_norm) && res_try_norm < last_residual
                aux.b .= aux_try.b
                aux.b_bar .= aux_try.b_bar
                aux.c .= aux_try.c
                aux.c_bar .= aux_try.c_bar
                accepted = true
                break
            end
        end
        if !accepted
            return JKSSolution(aux, iter, last_residual, false)
        end
    end
    return JKSSolution(aux, maxiter, last_residual, false)
end

end  # module Hubbard1DJKSNLIE
