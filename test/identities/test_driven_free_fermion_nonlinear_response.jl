# ─────────────────────────────────────────────────────────────────────────────
# Test: exact nonlinear response of the ac-driven free-fermion tight-binding
# chain, checked against an INDEPENDENT real-time simulation.
#
# The closed form under test (src: TightBindingV1D_driven.jl) is the Jacobi–Anger
# Bessel spectrum of the Peierls-driven single-band current
#
#     j(k, τ) = 2t sin(k + K sin ωτ),        K = E₀/ω  (dimensionless drive),
#
#   dc / 0-th:      2t J₀(K) sin k                            (dynamic localization)
#   n-th harmonic:  4t Jₙ(K) · (sin k if n even, cos k if n odd).
#
# The INDEPENDENT route builds the driven chain's current time series FROM
# SCRATCH — an RK4 integration of the time-dependent Schrödinger / von-Neumann
# equation with the Peierls Hamiltonian — and a bare discrete Fourier transform
# extracts each harmonic amplitude.  NOTHING in that route knows about Bessel
# functions or the Jacobi–Anger identity; if the DFT of the simulated current
# reproduces `driven_band_harmonic_weights` (= Bessel Jₙ), the closed form is
# confirmed by genuinely independent means.
#
# Pillars:
#   A. single-mode harmonic spectrum:  DFT of RK4 current == Bessel weights, for
#      k = π/2 (dc + even harmonics) and k = 0 (odd harmonics) — together these
#      pin every Jₙ(K).
#   B. dynamic localization (Dunlap–Kenkre):  cycle-averaged mode current
#      == 2t J₀(K) sin k across many k, and the whole band localizes at once at
#      the first Bessel zero K = 2.404826; ties to fetch(_, DynamicLocalization).
#   C. perturbative order-n limit:  small K ⇒ n-th harmonic ∝ Kⁿ (χ⁽ⁿ⁾), n = 1
#      recovers the linear conductivity slope 2t.
#   D. many-body integrator:  RK4 of the half-filled correlation matrix
#      dC/dτ = -i[H(τ),C] reproduces the summed single-mode current — exercises
#      the L×L matrix path and the additivity of the nonlinear response.
#   E. domain guards.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra
using QAtlas:
    TightBindingV1D, DynamicLocalization, Infinite, fetch, driven_band_harmonic_weights

# ── independent machinery: driven ring, RK4, bare DFT ────────────────────────
# (underscore-prefixed to stay clear of the shared included-test namespace)

# Vector potential of E(τ) = E₀ cos ωτ:  A(τ) = -(E₀/ω) sin ωτ = -K sin ωτ.
_dffnlr_A(τ, K, ω) = -K * sin(ω * τ)

# Peierls hopping Hamiltonian on an L-site ring: H[j,j+1] = -t e^{-iA}, h.c.
function _dffnlr_H(L, t, A)
    H = zeros(ComplexF64, L, L)
    ph = cis(-A)                       # e^{-iA}
    @inbounds for j in 1:L
        k = mod1(j + 1, L)
        H[j, k] += -t * ph
        H[k, j] += -t * conj(ph)
    end
    return H
end

# Uniform current operator J = -∂H/∂A: J[j,j+1] = -i t e^{-iA}, J[j+1,j] = +i t e^{+iA}.
function _dffnlr_J(L, t, A)
    J = zeros(ComplexF64, L, L)
    ph = cis(-A)
    @inbounds for j in 1:L
        k = mod1(j + 1, L)
        J[j, k] += -im * t * ph
        J[k, j] += im * t * conj(ph)
    end
    return J
end

# RK4 for i∂τψ = H(τ)ψ  ⟹  ∂τψ = -i H(τ)ψ.
function _dffnlr_step_psi(ψ, τ, dt, L, t, K, ω)
    f(τ, ψ) = -im * (_dffnlr_H(L, t, _dffnlr_A(τ, K, ω)) * ψ)
    k1 = f(τ, ψ)
    k2 = f(τ + dt / 2, ψ .+ (dt / 2) .* k1)
    k3 = f(τ + dt / 2, ψ .+ (dt / 2) .* k2)
    k4 = f(τ + dt, ψ .+ dt .* k3)
    return ψ .+ (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

# Current of a single Bloch mode |k0⟩ over one drive period: N samples at τ_m = m T/N.
function _dffnlr_series(L, t, K, ω, k0, N, S)
    T = 2π / ω
    dt = T / (N * S)
    ψ = ComplexF64[cis(k0 * j) for j in 1:L] ./ sqrt(L)
    τ = 0.0
    j = zeros(Float64, N)
    for m in 0:(N - 1)
        Jτ = _dffnlr_J(L, t, _dffnlr_A(τ, K, ω))
        j[m + 1] = real(ψ' * Jτ * ψ)
        for _ in 1:S
            ψ = _dffnlr_step_psi(ψ, τ, dt, L, t, K, ω)
            τ += dt
        end
    end
    return j
end

# Bare DFT harmonic amplitude of a period-sampled real signal.
function _dffnlr_amp(sig, n)
    N = length(sig)
    n == 0 && return sum(sig) / N                       # signed dc component
    c = 2 / N * sum(sig[m + 1] * cos(2π * n * m / N) for m in 0:(N - 1))
    s = 2 / N * sum(sig[m + 1] * sin(2π * n * m / N) for m in 0:(N - 1))
    return hypot(c, s)                                  # harmonic magnitude
end

# RK4 of the half-filled correlation matrix dC/dτ = -i[H(τ),C]; total current tr(J C).
function _dffnlr_series_corr(L, t, K, ω, C0, N, S)
    T = 2π / ω
    dt = T / (N * S)
    C = ComplexF64.(C0)
    fC(τ, C) = (H=_dffnlr_H(L, t, _dffnlr_A(τ, K, ω)); -im .* (H * C .- C * H))
    τ = 0.0
    j = zeros(Float64, N)
    for m in 0:(N - 1)
        Jτ = _dffnlr_J(L, t, _dffnlr_A(τ, K, ω))
        j[m + 1] = real(tr(Jτ * C))
        for _ in 1:S
            k1 = fC(τ, C)
            k2 = fC(τ + dt / 2, C .+ (dt / 2) .* k1)
            k3 = fC(τ + dt / 2, C .+ (dt / 2) .* k2)
            k4 = fC(τ + dt, C .+ dt .* k3)
            C = C .+ (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
            τ += dt
        end
    end
    return j
end

const _DFFNLR_J0_ZERO1 = 2.404825557695773   # first zero of J₀ (dynamic localization)

@testset "driven free-fermion nonlinear response — exact Bessel harmonics" begin
    t = 1.0
    ω = 1.0
    L = 16                       # k = 0 (n=0) and k = π/2 (n=L/4) are ring momenta
    N = 128                      # samples per drive period (Nyquist ≫ nmax)
    S = 8                        # RK4 substeps between samples

    # ── Pillar A — single-mode harmonic spectrum == Bessel weights ───────────
    @testset "A: harmonic spectrum == Bessel Jₙ(K)  (K=$K)" for K in (1.0, 2.0)
        w = driven_band_harmonic_weights(K; nmax=5)          # w[n+1] = Jₙ(K)

        # k = π/2:  sin k = 1, cos k = 0  ⇒  dc + even harmonics, weight 2t / 4t.
        jhalf = _dffnlr_series(L, t, K, ω, π / 2, N, S)
        @test isapprox(_dffnlr_amp(jhalf, 0), 2t * w[1]; atol=2e-3)      # dc = 2t J₀
        @test isapprox(_dffnlr_amp(jhalf, 2), 4t * abs(w[3]); atol=2e-3) # 2ω = 4t J₂
        @test isapprox(_dffnlr_amp(jhalf, 4), 4t * abs(w[5]); atol=2e-3) # 4ω = 4t J₄
        @test _dffnlr_amp(jhalf, 1) < 2e-3                               # odd absent
        @test _dffnlr_amp(jhalf, 3) < 2e-3

        # k = 0:  sin k = 0, cos k = 1  ⇒  odd harmonics only, weight 4t.
        jzero = _dffnlr_series(L, t, K, ω, 0.0, N, S)
        @test isapprox(_dffnlr_amp(jzero, 1), 4t * abs(w[2]); atol=2e-3) # 1ω = 4t J₁
        @test isapprox(_dffnlr_amp(jzero, 3), 4t * abs(w[4]); atol=2e-3) # 3ω = 4t J₃
        @test isapprox(_dffnlr_amp(jzero, 5), 4t * abs(w[6]); atol=2e-3) # 5ω = 4t J₅
        @test abs(_dffnlr_amp(jzero, 0)) < 2e-3                          # dc absent
        @test _dffnlr_amp(jzero, 2) < 2e-3                               # even absent
    end

    # ── Pillar B — dynamic localization (Dunlap–Kenkre 1986) ─────────────────
    @testset "B: t_eff = t·J₀(K) and band collapse at first Bessel zero" begin
        m = TightBindingV1D(; t=t)                      # free-fermion point (V=0)

        # closed-form registered quantity: t_eff = t J₀(K)
        w0(K) = driven_band_harmonic_weights(K; nmax=0)[1]
        for K in (0.0, 0.8, 1.5, 2.0, _DFFNLR_J0_ZERO1)
            @test isapprox(
                fetch(m, DynamicLocalization(), Infinite(); drive=K), t * w0(K); atol=1e-12
            )
        end
        @test fetch(m, DynamicLocalization(), Infinite(); drive=0.0) == t   # no drive
        @test abs(fetch(m, DynamicLocalization(), Infinite(); drive=_DFFNLR_J0_ZERO1)) <
            1e-6                     # localized

        # independent: cycle-averaged mode current == 2t J₀(K) sin k, all modes.
        for K in (1.0, 2.0), n in (1, 3, 5, 7)          # k = 2π n / L
            k0 = 2π * n / L
            dc = _dffnlr_amp(_dffnlr_series(L, t, K, ω, k0, N, S), 0)
            @test isapprox(dc, 2t * w0(K) * sin(k0); atol=3e-3)
        end
        # at K = 2.4048 every mode's dc current collapses → whole band localizes.
        for n in 1:7
            k0 = 2π * n / L
            dc = _dffnlr_amp(_dffnlr_series(L, t, _DFFNLR_J0_ZERO1, ω, k0, N, S), 0)
            @test abs(dc) < 3e-3
        end
    end

    # ── Pillar C — small drive: n-th harmonic is the order-n response ─────────
    @testset "C: perturbative χ⁽ⁿ⁾ scaling and linear-response slope" begin
        Ksmall = 0.15
        w = driven_band_harmonic_weights(Ksmall; nmax=3)
        # Jₙ(K) ≈ (K/2)ⁿ / n!  ⇒  successive harmonics shrink by ~ K/2.
        @test isapprox(w[2], Ksmall / 2; rtol=2e-2)              # J₁ ≈ K/2
        @test isapprox(w[3], (Ksmall / 2)^2 / 2; rtol=3e-2)      # J₂ ≈ (K/2)²/2
        # linear-response slope of the first harmonic (k=0): amp₁ = 4t J₁ ≈ 2t·K.
        amp1 = _dffnlr_amp(_dffnlr_series(L, t, Ksmall, ω, 0.0, N, S), 1)
        @test isapprox(amp1 / Ksmall, 2t; rtol=2e-2)
        # order separation: 2nd harmonic is O(K) smaller than the 1st.
        amp2 = _dffnlr_amp(_dffnlr_series(L, t, Ksmall, ω, π / 2, N, S), 2)
        @test amp2 < amp1 * Ksmall
    end

    # ── Pillar D — many-body correlation-matrix integrator ───────────────────
    @testset "D: half-filled correlation-matrix RK4 == summed single-mode" begin
        Ld = 12
        K = 1.3
        Nd, Sd = 96, 8
        kk = [2π * n / Ld for n in 0:(Ld - 1)]
        occ = partialsortperm([-2t * cos(k) for k in kk], 1:(Ld ÷ 2))  # lowest ½
        kocc = kk[occ]
        # C₀ = Σ_{k∈occ} |k⟩⟨k|,  |k⟩ = L^{-1/2} e^{ikj}
        C0 = zeros(ComplexF64, Ld, Ld)
        for k in kocc, a in 1:Ld, b in 1:Ld
            C0[a, b] += cis(k * (a - b)) / Ld
        end
        jrk4 = _dffnlr_series_corr(Ld, t, K, ω, C0, Nd, Sd)
        T = 2π / ω
        jref = [
            sum(2t * sin(k - _dffnlr_A(m * T / Nd, K, ω)) for k in kocc) for m in 0:(Nd - 1)
        ]
        @test maximum(abs, jrk4 .- jref) < 5e-3
    end

    # ── Pillar E — domain guards ─────────────────────────────────────────────
    @testset "E: guards" begin
        @test_throws DomainError fetch(
            TightBindingV1D(; t=1.0, V=0.5), DynamicLocalization(), Infinite(); drive=1.0
        )
        @test_throws DomainError driven_band_harmonic_weights(1.0; nmax=-1)
        @test length(driven_band_harmonic_weights(1.0)) == 7    # nmax=6 default ⇒ 0..6
    end
end
