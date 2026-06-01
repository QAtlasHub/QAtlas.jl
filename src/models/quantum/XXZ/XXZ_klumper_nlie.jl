# ─────────────────────────────────────────────────────────────────────────────
# src/models/quantum/XXZ/XXZ_klumper_nlie.jl
#
# Klümper Quantum-Transfer-Matrix (QTM) non-linear integral equation (NLIE)
# for the critical XXZ chain (-1 < Δ < 1) in the thermodynamic limit, at
# arbitrary temperature and zero magnetic field.
#
# Following Klümper, Z. Phys. B 91, 507 (1993)  [DOI 10.1007/BF01316831],
# equations (5.4)-(5.7).  The Hamiltonian convention is matched to
# QAtlas's XXZ1D
#
#   H = J Σᵢ [ Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Δ Sᶻᵢ Sᶻᵢ₊₁ ]
#
# via γ = arccos(Δ), β̃ = β J sin(γ)/2 (pinned by the XX (Δ=0) limit
# against the exact free-fermion result; see _xx_limit_check in the
# test file).  In Klümper's energy units the free-energy density per
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
function _kshift(x::Real, γ::Real, ε::Real;
                 kmax::Real=max(60.0, 30.0 / ε),
                 rtol::Real=1e-10)
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

function build_grid(γ::Real;
                    N::Int=512,
                    L_factor::Real=20.0,
                    ε_shift::Real=0.5,
                    rtol::Real=1e-10)
    0 < γ < π || throw(DomainError(γ, "γ must lie in (0, π) for the critical regime"))
    L = L_factor * γ / π
    dx = 2L / N
    x = [(i - (N + 1) / 2) * dx for i in 1:N]
    sech_factor = [1 / cosh(π * xi / γ) for xi in x]
    driving_template = [-π * s / γ for s in sech_factor]
    Nk = 2N - 1
    k_vec = Vector{Float64}(undef, Nk)
    k_shift_vec = Vector{ComplexF64}(undef, Nk)
    for (i, d) in enumerate(-(N - 1):(N - 1))
        xd = d * dx
        k_vec[i] = _kreal(xd, γ; rtol=rtol)
        k_shift_vec[i] = _kshift(xd, γ, ε_shift; rtol=rtol)
    end
    return NLIEGrid(Float64(γ), N, dx, L, Float64(ε_shift),
                    x, sech_factor, driving_template, k_vec, k_shift_vec)
end

# Toeplitz matrix-vector product g[i] = Σⱼ k_vec[N + i - j] · f[j] · dx.
function _conv!(g::AbstractVector, k_vec::AbstractVector,
                f::AbstractVector, N::Int, dx::Real)
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

function solve_klumper_nlie(grid::NLIEGrid, β̃::Real;
                            α::Real=0.5,
                            maxiter::Int=500,
                            tol::Real=1e-10,
                            verbose::Bool=false)
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
    return γ, β̃, escale
end

end  # module XXZKlumperNLIE

# ════════════════════════════════════════════════════════════════════════════
# Cache + dispatch glue (outside the module so it can call into QAtlas)
# ════════════════════════════════════════════════════════════════════════════

const _XXZ_NLIE_GRID_CACHE = Dict{Tuple{Float64,Int,Float64,Float64},XXZKlumperNLIE.NLIEGrid}()

function _xxz_nlie_grid(γ::Float64; N::Int=128, L_factor::Float64=15.0, ε_shift::Float64=0.5)
    key = (round(γ; digits=10), N, L_factor, ε_shift)
    haskey(_XXZ_NLIE_GRID_CACHE, key) && return _XXZ_NLIE_GRID_CACHE[key]
    g = XXZKlumperNLIE.build_grid(γ; N=N, L_factor=L_factor, ε_shift=ε_shift)
    _XXZ_NLIE_GRID_CACHE[key] = g
    return g
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
    if !(-0.99 < Δ < 0.99)
        @warn "XXZ1D Klümper NLIE skips |Δ| ≥ 0.99 (mapping degenerates as sin γ → 0); " *
            "Heisenberg/FM endpoint deferred to issue #521." Δ
        return NaN
    end
    γ, β̃, escale = XXZKlumperNLIE.qatlas_to_klumper(Float64(Δ), Float64(J), Float64(beta))
    grid = _xxz_nlie_grid(γ)
    sol = XXZKlumperNLIE.solve_klumper_nlie(grid, β̃; α=0.4, maxiter=400, tol=1e-9)
    if !sol.converged
        @warn "XXZ1D Klümper NLIE did not converge" Δ beta residual=sol.residual iterations=sol.iterations
        return NaN
    end
    return escale * XXZKlumperNLIE.free_energy_excess_klumper(sol)
end
