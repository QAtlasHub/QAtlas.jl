# Verification of the AKLT1D biquadratic-aware finite-T HTSE (#506).
#
# Independent references (never the src formulas re-typed):
#   • per-bond cumulants κₙ recomputed from exact ED moments Tr(Hᵏ)/3ᴺ;
#   • the bilinear-limit anchor κ₂ = r²/3 = 4/3 — the published per-bond
#     specific-heat coefficient d₂ of Lohmann–Schmidt–Richter 2014 (r=s(s+1)=2);
#   • a full Boltzmann-ED of a finite PBC ring (Z = Σ e^{-βE}, a wholly
#     different computation from the cumulant route);
#   • autodiff for the c_v = β²φ'' and s = β(ε−f) thermodynamic identities;
#   • the β→0 limits (s→ln3, c→0, f→−ln3/β).

using QAtlas, Test, LinearAlgebra
using QAtlas:
    AKLT1D, Infinite, FreeEnergy, SpecificHeat, ThermalEntropy, SusceptibilityZZ, fetch
using ForwardDiff: ForwardDiff

# ── local spin-1 bilinear-biquadratic builders (independent reference) ──────
const _Sz = Matrix{ComplexF64}([1.0 0 0; 0 0 0; 0 0 -1])
const _Sp = ComplexF64[0 sqrt(2) 0; 0 0 sqrt(2); 0 0 0]
const _Sx = (_Sp + _Sp') / 2
const _Sy = (_Sp - _Sp') / (2im)
const _SdotS = real(kron(_Sx, _Sx) + kron(_Sy, _Sy) + kron(_Sz, _Sz))
_idn(n) = Matrix{Float64}(I, 3^n, 3^n)

# OBC chain H = Σ J1 S·S + J2 (S·S)²
function _blbq_obc(N, J1, J2)
    H = zeros(Float64, 3^N, 3^N)
    bond = J1 * _SdotS + J2 * _SdotS^2
    for i in 1:(N - 1)
        H += kron(_idn(i - 1), kron(bond, _idn(N - i - 1)))
    end
    return H
end

# add the periodic wrap bond N--1
function _blbq_pbc(N, J1, J2)
    H = _blbq_obc(N, J1, J2)
    Sw = real(sum(kron(Sa, kron(_idn(N - 2), Sa)) for Sa in (_Sx, _Sy, _Sz)))
    return H + J1 * Sw + J2 * Sw^2
end

# extensive infinite-T cumulants from raw moments mₖ = Tr(Hᵏ)/3ᴺ
function _cums(H)
    d = size(H, 1)
    m = [tr(H^k) / d for k in 1:4]
    return [
        m[1],
        m[2] - m[1]^2,
        m[3] - 3m[1] * m[2] + 2m[1]^3,
        m[4] - 4m[1] * m[3] - 3m[2]^2 + 12m[1]^2 * m[2] - 6m[1]^4,
    ]
end

# per-bond cumulant = thermodynamic-limit increment Cₙ(N) − Cₙ(N−1)
_perbond_kappa(J2; N=6) = _cums(_blbq_obc(N, 1.0, J2)) .- _cums(_blbq_obc(N - 1, 1.0, J2))

# independent Boltzmann-ED per-site thermodynamics
function _ed_thermo(E, β, N)
    Emin = minimum(E)
    w = exp.(-β .* (E .- Emin))
    Z = sum(w)
    Ea = sum(E .* w) / Z
    E2 = sum((E .^ 2) .* w) / Z
    φ = (log(Z) - β * Emin) / N
    return (cv=β^2 * (E2 - Ea^2) / N, f=(-φ / β), s=φ + β * Ea / N)
end

@testset "AKLT1D biquadratic-aware HTSE (#506)" begin
    @testset "hardcoded κ vs exact ED moments (converged in N)" begin
        κ_ed = _perbond_kappa(1 / 3; N=6)
        @test κ_ed ≈ collect(QAtlas._AKLT_HTSE_KAPPA) atol = 1e-9
        @test _perbond_kappa(1 / 3; N=5) ≈ κ_ed atol = 1e-9   # converged
    end

    @testset "bilinear limit κ₂ = r²/3 = 4/3  [Lohmann2014 per-bond d₂]" begin
        κ_bil = _perbond_kappa(0.0)
        @test κ_bil[1] ≈ 0.0 atol = 1e-10        # Tr(S·S) = 0
        @test κ_bil[2] ≈ 4 / 3 atol = 1e-10      # r = s(s+1) = 2 → r²/3
    end

    @testset "HTSE vs independent Boltzmann-ED (N=7 PBC), βJ ≲ 0.35" begin
        m = AKLT1D(; J=1.0)
        E = eigvals(Symmetric(_blbq_pbc(7, 1.0, 1 / 3)))
        for β in (0.1, 0.2, 0.3, 0.35)
            ed = _ed_thermo(E, β, 7)
            @test fetch(m, SpecificHeat(), Infinite(); scheme=:htse, beta=β) ≈ ed.cv rtol =
                0.015
            @test fetch(m, FreeEnergy(), Infinite(); scheme=:htse, beta=β) ≈ ed.f rtol =
                2e-3
            @test fetch(m, ThermalEntropy(), Infinite(); scheme=:htse, beta=β) ≈ ed.s rtol =
                4e-3
        end
    end

    @testset "thermodynamic identities via autodiff (c_v=β²φ'', s=β(ε−f))" begin
        m = AKLT1D(; J=0.8)
        φ(β) = -β * fetch(m, FreeEnergy(), Infinite(); scheme=:htse, beta=β)  # lnZ/N
        for β in (0.15, 0.3)
            cv_ad = β^2 * ForwardDiff.derivative(b -> ForwardDiff.derivative(φ, b), β)
            @test fetch(m, SpecificHeat(), Infinite(); scheme=:htse, beta=β) ≈ cv_ad rtol =
                1e-8
            ε = -ForwardDiff.derivative(φ, β)
            f = fetch(m, FreeEnergy(), Infinite(); scheme=:htse, beta=β)
            s = fetch(m, ThermalEntropy(), Infinite(); scheme=:htse, beta=β)
            @test s ≈ β * (ε - f) rtol = 1e-9
        end
    end

    @testset "high-T limits (β → 0)" begin
        m = AKLT1D(; J=1.0)
        @test fetch(m, ThermalEntropy(), Infinite(); scheme=:htse, beta=1e-4) ≈ log(3) atol =
            1e-6
        @test fetch(m, SpecificHeat(), Infinite(); scheme=:htse, beta=1e-4) < 1e-6
        @test fetch(m, FreeEnergy(), Infinite(); scheme=:htse, beta=1e-3) ≈ -log(3) / 1e-3 rtol =
            1e-3
    end

    @testset "J-scaling: c_v(J, β) = c_v(1, Jβ)" begin
        @test fetch(AKLT1D(; J=2.0), SpecificHeat(), Infinite(); scheme=:htse, beta=0.15) ≈
            fetch(AKLT1D(; J=1.0), SpecificHeat(), Infinite(); scheme=:htse, beta=0.30) rtol =
            1e-12
    end

    @testset "scheme routing + domain guards" begin
        m = AKLT1D(; J=1.0)
        # canonical stays the exact β=∞ limit
        @test fetch(m, FreeEnergy(), Infinite(); beta=Inf) ≈ -2 / 3
        @test fetch(m, SpecificHeat(), Infinite(); scheme=:canonical, beta=Inf) == 0.0
        # :htse rejects β=∞ and β≤0
        @test_throws DomainError fetch(m, FreeEnergy(), Infinite(); scheme=:htse, beta=Inf)
        @test_throws DomainError fetch(
            m, SpecificHeat(), Infinite(); scheme=:htse, beta=-1.0
        )
        # unknown scheme → ArgumentError; finite-β χ HTSE not implemented (#506 follow-up)
        @test_throws ArgumentError fetch(
            m, FreeEnergy(), Infinite(); scheme=:bogus, beta=1.0
        )
        @test_throws DomainError fetch(m, SusceptibilityZZ(), Infinite(); beta=1.0)
    end
end
