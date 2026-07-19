# ─────────────────────────────────────────────────────────────────────────────
# src/models/quantum/XXZ/XXZ_klumper_nlie.jl
#
# Klümper Quantum-Transfer-Matrix (QTM) non-linear integral equation (NLIE)
# for the critical XXZ chain (-1 < Δ < 1) in the thermodynamic limit, at
# arbitrary temperature and zero magnetic field.
#
# Following Klümper, [Klumper1993](@cite)  [DOI 10.1007/BF01316831],
# equations (5.4)-(5.7).  The Hamiltonian convention is matched to
# QAtlas's XXZ1D
#
#   H = J Σᵢ [ Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Δ Sᶻᵢ Sᶻᵢ₊₁ ]
#
# via γ = arccos(Δ), β̃ = β J sin(γ)/2.  The factor sin(γ)/2 is the
# rescaling that maps Klümpers natural energy unit (eq. 5.1, J_X = 1/sin γ
# after the canonical sign-flip and axis-cycling unitaries) onto QAtlass
# H_Q = J Σ S·S coupling; the prefactor sin γ is recovered explicitly in
# the XX limit (γ = π/2, eq. 5.7), where the spectral parameter identifi-
# cation reduces to β̃ = β J / 2 and is cross-checked against the
# free-fermion closed form by the XX-limit testset.  In Klümper's energy units the free-energy density per
# site is
#
#       β̃ ( f_K − ε_K ) = − (1/(2γ)) ∫_{-∞}^{∞} ln[ A(x) Ā(x) ] / cosh(πx/γ) dx.   (5.5)
#
# Converting to QAtlas units multiplies by (J sin(γ)/2):
#
#       f_Q(T) − ε_0_Q = (J sin(γ)/2)( f_K − ε_K ).
#
# Validity. The mapping degenerates at sin(γ) → 0 (Δ → ±1: the
# Heisenberg and FM-Heisenberg points). We refuse those endpoints
# at the dispatch boundary and leave them as follow-up (issue #521).
#
# Numerical strategy. Picard iteration with mixing α on a uniform
# grid in x ∈ [-L, L]. Convolutions are evaluated as direct Toeplitz
# matrix-vector products (O(N²) per iteration; FFTW intentionally
# avoided to keep the dependency surface minimal). The shifted
# kernel k(x - iγ + iε) is regularised with a small ε > 0 so its
# Fourier coefficient ŝ(k) e^{(γ-ε)|k|} on k > 0 decays as e^{-εk}.
# ─────────────────────────────────────────────────────────────────────────────

"""
    XXZKlumperNLIE

Internal numerical solver for the Klümper Quantum-Transfer-Matrix
non-linear integral equation (NLIE) of the critical XXZ chain
(-1 < Δ < 1) in the thermodynamic limit.  Used by
`fetch(::XXZ1D, ::FreeEnergy, ::Infinite; beta=...)` to provide an
analytic reference free-energy density at general Δ beyond the
XX (Δ = 0) free-fermion point.

This is an implementation detail of the XXZ dispatch surface, not a
public API.  See the module-level header comment for the full equation
references (Klümper 1993, eq. 5.4-5.7).
"""
module XXZKlumperNLIE

using LinearAlgebra
using QuadGK

# ════════════════════════════════════════════════════════════════════════════
# Spectral kernel and its real-axis / shifted samples
# ════════════════════════════════════════════════════════════════════════════

# Klümper 1993, eq. (5.6): the Fourier-transformed kernel of the NLIE.
@inline function _khat(k::Real, γ::Real)
    iszero(k) && return zero(float(k))
    num = sinh((π / 2 - γ) * k)
    den = 2 * cosh(γ * k / 2) * sinh((π - γ) * k / 2)
    return num / den
end

# Real-axis sample k(x) = (1/π) ∫_0^∞ ŝ(k) cos(kx) dk (even function of x).
function _kreal(x::Real, γ::Real; kmax::Real=200.0, rtol::Real=1e-10)
    val, _ = quadgk(k -> _khat(k, γ) * cos(k * x), 0.0, kmax; rtol=rtol)
    return val / π
end

# Shifted kernel k(x - i(γ - ε)), complex-valued. Decay at large |k| is
# controlled by ε > 0 (must satisfy 0 < ε < γ).
# Adaptive kmax for the slow e^{-εk} tail of the iε-regularised kernel.
_kshift_kmax(ε::Real) = max(60.0, 30.0 / ε)

# Minimum safe ε_shift: keeps cosh(k·(γ-ε)) within Float64 range during
# numerical integration. Below this floor _kshift overflows quadgk.
const _KSHIFT_EPS_MIN = 0.1

function _kshift(x::Real, γ::Real, ε::Real; kmax::Real=_kshift_kmax(ε), rtol::Real=1e-10)
    decay = γ - ε
    decay > 0 || throw(ArgumentError("ε must satisfy 0 < ε < γ; got ε=$ε, γ=$γ"))
    re, _ = quadgk(k -> _khat(k, γ) * cos(k * x) * cosh(k * decay), 0.0, kmax; rtol=rtol)
    im, _ = quadgk(k -> _khat(k, γ) * sin(k * x) * sinh(k * decay), 0.0, kmax; rtol=rtol)
    return complex(re / π, im / π)
end

# ════════════════════════════════════════════════════════════════════════════
# Grid + precomputed Toeplitz kernel tables
# ════════════════════════════════════════════════════════════════════════════

struct NLIEGrid
    γ::Float64
    N::Int
    dx::Float64
    L::Float64
    ε_shift::Float64
    x::Vector{Float64}
    sech_factor::Vector{Float64}             # 1 / cosh(πx/γ)
    driving_template::Vector{Float64}        # −π / (γ cosh(πx/γ)), without β̃
    k_vec::Vector{Float64}                   # length 2N-1, k(d·dx)
    k_shift_vec::Vector{ComplexF64}          # length 2N-1, k_shift(d·dx)
end

function build_grid(
    γ::Real; N::Int=512, L_factor::Real=20.0, ε_shift::Real=0.5, rtol::Real=1e-10
)
    0 < γ < π || throw(DomainError(γ, "γ must lie in (0, π) for the critical regime"))
    ε_shift >= _KSHIFT_EPS_MIN || throw(
        ArgumentError(
            "ε_shift = $(ε_shift) below safe floor $(_KSHIFT_EPS_MIN); " *
            "the shifted kernel integral overflows for smaller ε. " *
            "Choose ε_shift in [$(_KSHIFT_EPS_MIN), γ) or use the default 0.5.",
        ),
    )
    ε_shift < γ ||
        throw(ArgumentError("ε_shift = $(ε_shift) must be strictly less than γ = $(γ)."))
    # γ-aware overflow guard. `_kshift_kmax(ε)` returns max(60, 30/ε); the
    # cosh(k·(γ-ε)) integrand evaluated at k = kmax overflows Float64 above
    # cosh(709). If kmax·(γ-ε) > 709, quadgk silently returns NaN and poisons
    # the entire Toeplitz kernel. Refuse such inputs with an explicit error.
    decay_c = γ - ε_shift
    kmax_guess = _kshift_kmax(ε_shift)
    if kmax_guess * decay_c > 700.0
        eps_floor_for_γ = γ - 700.0 / kmax_guess
        throw(
            ArgumentError(
                "ε_shift = $(ε_shift) combined with γ = $(γ) (Δ = $(cos(γ))) drives " *
                "cosh(kmax·(γ-ε)) over the Float64 ceiling; the shifted kernel " *
                "integral would silently return NaN. Increase ε_shift above " *
                "$(eps_floor_for_γ) or reduce γ.",
            ),
        )
    end
    L = L_factor * γ / π
    dx = 2L / N
    x = [(i - (N + 1) / 2) * dx for i in 1:N]
    sech_factor = [1 / cosh(π * xi / γ) for xi in x]
    driving_template = [-π * s / γ for s in sech_factor]
    Nk = 2N - 1
    k_vec = Vector{Float64}(undef, Nk)
    k_shift_vec = Vector{ComplexF64}(undef, Nk)
    for (i, d) in enumerate((-(N - 1)):(N - 1))
        xd = d * dx
        k_vec[i] = _kreal(xd, γ; rtol=rtol)
        k_shift_vec[i] = _kshift(xd, γ, ε_shift; rtol=rtol)
    end
    return NLIEGrid(
        Float64(γ),
        N,
        dx,
        L,
        Float64(ε_shift),
        x,
        sech_factor,
        driving_template,
        k_vec,
        k_shift_vec,
    )
end

# Toeplitz matrix-vector product g[i] = Σⱼ k_vec[N + i - j] · f[j] · dx.
function _conv!(
    g::AbstractVector, k_vec::AbstractVector, f::AbstractVector, N::Int, dx::Real
)
    @inbounds for i in 1:N
        s = zero(eltype(g))
        for j in 1:N
            s += k_vec[N + i - j] * f[j]
        end
        g[i] = s * dx
    end
    return g
end

# ════════════════════════════════════════════════════════════════════════════
# Solver (Picard with mixing)
# ════════════════════════════════════════════════════════════════════════════

struct NLIESolution
    grid::NLIEGrid
    β̃::Float64
    a::Vector{ComplexF64}
    converged::Bool
    iterations::Int
    residual::Float64
end

function solve_klumper_nlie(
    grid::NLIEGrid,
    β̃::Real;
    α::Real=0.5,
    maxiter::Int=500,
    tol::Real=1e-10,
    verbose::Bool=false,
)
    N = grid.N
    driving = β̃ .* grid.driving_template
    a = ComplexF64.(exp.(driving))
    ln_a = log.(a)
    ln_A = log.(1 .+ a)
    int1 = zeros(ComplexF64, N)
    int2 = zeros(ComplexF64, N)
    ln_a_new = zeros(ComplexF64, N)
    converged = false
    residual = Inf
    iter = 0
    for it in 1:maxiter
        _conv!(int1, grid.k_vec, ln_A, N, grid.dx)
        _conv!(int2, grid.k_shift_vec, conj.(ln_A), N, grid.dx)
        @. ln_a_new = driving + int1 - int2
        residual = maximum(abs.(ln_a_new .- ln_a))
        @. ln_a = (1 - α) * ln_a + α * ln_a_new
        a .= exp.(ln_a)
        ln_A .= log.(1 .+ a)
        iter = it
        verbose && it % 20 == 0 && @info "NLIE iter=$it res=$residual"
        if residual < tol
            converged = true
            break
        end
    end
    return NLIESolution(grid, Float64(β̃), copy(a), converged, iter, residual)
end

# Klümper-units free-energy excess: β̃ (f_K − ε_K), eq. (5.5).
function free_energy_excess_klumper(sol::NLIESolution)
    A = 1 .+ sol.a
    Abar = 1 .+ conj.(sol.a)
    integrand = real.(log.(A) .+ log.(Abar)) .* sol.grid.sech_factor
    integral = sum(integrand) * sol.grid.dx
    return -integral / (2 * sol.grid.γ * sol.β̃)
end

# ════════════════════════════════════════════════════════════════════════════
# QAtlas convention map
# ════════════════════════════════════════════════════════════════════════════

# Returns (γ, β̃, escale) where escale converts Klümper free-energy
# excess to QAtlas units: (f_Q − ε_0_Q) = escale · (f_K − ε_K).
function qatlas_to_klumper(Δ::Real, J::Real, β::Real)
    γ = acos(Δ)
    sinγ = sin(γ)
    β̃ = β * J * sinγ / 2
    escale = J * sinγ / 2
    return (; γ, β̃, escale)
end

end  # module XXZKlumperNLIE

# ════════════════════════════════════════════════════════════════════════════
# Cache + dispatch glue (outside the module so it can call into QAtlas)
# ════════════════════════════════════════════════════════════════════════════

# Module-level grid cache. Practical bound: one entry per distinct
# (γ, N, L_factor, ε_shift) tuple. In production usage all calls go through
# `_xxz_nlie_grid`, which fixes (N, L_factor, ε_shift) = (128, 15.0, 0.5),
# so the cache size grows by one per distinct Δ (= arccos γ, rounded to 10
# digits). Memory per entry ≈ 4·N + 2·(2N-1)·8 bytes ≈ 4 KiB; a 100-Δ
# coverage sweep stays under 0.5 MiB.
const _XXZ_NLIE_GRID_CACHE = Dict{
    Tuple{Float64,Int,Float64,Float64},XXZKlumperNLIE.NLIEGrid
}()
const _XXZ_NLIE_GRID_CACHE_LOCK = ReentrantLock()

# `_xxz_nlie_grid` returns a cached `NLIEGrid` keyed on
# `(round(γ; digits=10), N, L_factor, ε_shift)`. N defaults to 128 —
# the production sweet spot empirically: at γ ≈ π/3 the kernel decays
# as e^{-π|x|/γ}, so L_factor·γ/π = 15·1/3 ≈ 5 with dx ≈ 0.04 resolves
# the kernel to ~1e-7. Coarser than the standalone `build_grid` default
# (N=512), but the dispatch-level caller cares about wall time over the
# spectral tail beyond 5γ/π, which is already e^{-15}-suppressed.
function _xxz_nlie_grid(
    γ::Float64; N::Int=128, L_factor::Float64=15.0, ε_shift::Float64=0.5
)
    key = (round(γ; digits=10), N, L_factor, ε_shift)
    lock(_XXZ_NLIE_GRID_CACHE_LOCK) do
        get!(_XXZ_NLIE_GRID_CACHE, key) do
            return XXZKlumperNLIE.build_grid(γ; N=N, L_factor=L_factor, ε_shift=ε_shift)
        end
    end
end

"""
    _xxz_klumper_free_energy_excess(model, beta) -> Float64

Per-site `f(β) − f(∞)` for `XXZ1D` in the critical regime (`-1 < Δ < 1`,
both endpoints excluded) at zero magnetic field, computed from the
Klümper QTM NLIE. Returns `NaN` with a warning at the endpoints or
on iteration failure.
"""
function _xxz_klumper_free_energy_excess(model::XXZ1D, beta::Real)
    Δ, J = model.Δ, model.J
    beta > 0 ||
        throw(DomainError(beta, "XXZ1D Klümper NLIE requires β > 0; got β = $(beta)."))
    J > 0 || throw(
        DomainError(
            J,
            "XXZ1D Klümper NLIE requires J > 0; got J = $(J). " *
            "J = 0 makes the Klümper energy unit J·sin(γ)/2 vanish (β̃ = 0), " *
            "causing a 0/0 in the free-energy integral.",
        ),
    )
    if !(-0.99 < Δ < 0.99)
        @warn "XXZ1D Klümper NLIE skips |Δ| ≥ 0.99 (mapping degenerates as sin γ → 0); " *
            "Heisenberg/FM endpoint deferred to issue #521." Δ
        return NaN
    end
    pars = XXZKlumperNLIE.qatlas_to_klumper(Float64(Δ), Float64(J), Float64(beta))
    grid = _xxz_nlie_grid(pars.γ)
    # Solver parameters are hard-coded here for the production dispatch path.
    # Callers needing tighter convergence (tol=1e-12) or different mixing α
    # should call XXZKlumperNLIE.solve_klumper_nlie(grid, β̃; α=..., ...) directly.
    sol = XXZKlumperNLIE.solve_klumper_nlie(grid, pars.β̃; α=0.4, maxiter=400, tol=1e-9)
    if !sol.converged
        @warn "XXZ1D Klümper NLIE did not converge for FreeEnergy@Infinite" Δ beta residual=sol.residual iterations=sol.iterations
        return NaN
    end
    # Surface near-marginal solutions even when the boolean flag is true.
    # `@info` rather than `@warn`: the result is mildly noisier than the
    # tol=1e-9 target, not wrong, so we avoid diluting the genuine
    # non-convergence warning above by reserving `@warn` for that.
    if sol.residual > 1e-7
        @info "XXZ1D Klümper NLIE converged but residual is above 1e-7; the result may be slightly biased." Δ beta residual=sol.residual
    end
    return pars.escale * XXZKlumperNLIE.free_energy_excess_klumper(sol)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Thermal entropy & specific heat via central finite difference of f(β)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Thermodynamic identities (all at zero magnetic field):
#
#   u(β) = ∂(β f) / ∂β = f(β) + β f'(β)
#   s(β) = β² f'(β)                            (entropy density)
#   c(β) = -β · ∂s/∂β = -2 β² f'(β) - β³ f''(β) (heat capacity)
#
# We compute f at 3 β values via the cached NLIE grid (one full grid build
# at the first call; subsequent calls at the same γ reuse the cache, so
# each extra β costs ~one Picard iteration, ~1–2 s in production).
#
# The relative step δrel = δ/β has a sweet spot near 1e-4: too small and
# round-off in f dominates the second difference; too large and the
# discretization error in f''(β) is no longer leading. We expose it as
# a kwarg so callers needing higher precision can tighten.

"""
    _xxz_klumper_thermal_derivatives(model, beta; δrel=1e-4) -> (Float64, Float64)

Return `(s, c)` from a 3-point central finite difference of
`_xxz_klumper_free_energy_excess`. NaN+warn for `|Δ| ≥ 0.99`, identical
to the FreeEnergy gate. Both `s` and `c` propagate NaN if any of the
three f evaluations fail to converge.
"""
function _xxz_klumper_thermal_derivatives(model::XXZ1D, beta::Real; δrel::Real=1e-4)
    Δ, J = model.Δ, model.J
    beta > 0 || throw(
        DomainError(
            beta, "XXZ1D Klümper thermal derivatives require β > 0; got β = $(beta)."
        ),
    )
    J > 0 || throw(
        DomainError(J, "XXZ1D Klümper thermal derivatives require J > 0; got J = $(J).")
    )
    δrel > 0 || throw(DomainError(δrel, "δrel must be > 0; got δrel = $(δrel)."))
    if !(-0.99 < Δ < 0.99)
        @warn "XXZ1D Klümper NLIE skips |Δ| ≥ 0.99 for thermal derivatives" Δ
        return (NaN, NaN)
    end
    δ = beta * δrel
    f_minus = _xxz_klumper_free_energy_excess(model, beta - δ)
    f_zero = _xxz_klumper_free_energy_excess(model, beta)
    f_plus = _xxz_klumper_free_energy_excess(model, beta + δ)
    if isnan(f_minus) || isnan(f_zero) || isnan(f_plus)
        return (NaN, NaN)
    end
    fp = (f_plus - f_minus) / (2δ)
    fpp = (f_plus - 2 * f_zero + f_minus) / δ^2
    s = beta^2 * fp
    c = -2 * beta^2 * fp - beta^3 * fpp
    return (s, c)
end

"""
    _xxz_klumper_entropy(model, beta; δrel=1e-4) -> Float64

Per-site Gibbs entropy density `s(β) = β² · ∂f/∂β` via central finite
difference of the Klümper NLIE FreeEnergy. NaN+warn at |Δ| ≥ 0.99.
"""
function _xxz_klumper_entropy(model::XXZ1D, beta::Real; δrel::Real=1e-4)
    return _xxz_klumper_thermal_derivatives(model, beta; δrel=δrel)[1]
end

"""
    _xxz_klumper_specific_heat(model, beta; δrel=1e-4) -> Float64

Per-site heat capacity `c(β) = -2 β² f'(β) - β³ f''(β)` via central
finite difference of the Klümper NLIE FreeEnergy. NaN+warn at |Δ| ≥ 0.99.
"""
function _xxz_klumper_specific_heat(model::XXZ1D, beta::Real; δrel::Real=1e-4)
    return _xxz_klumper_thermal_derivatives(model, beta; δrel=δrel)[2]
end
