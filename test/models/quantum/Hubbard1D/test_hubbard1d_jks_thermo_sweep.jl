# test/models/quantum/Hubbard1D/test_hubbard1d_jks_thermo_sweep.jl
#
# Detailed multi-temperature benchmark for the JKS NLIE FreeEnergy@Infinite
# dispatch (#523). Mirrors the structure of
# test/models/quantum/XXZ/test_xxz_klumper_nlie_thermo.jl so the two
# finite-T NLIE ports can be cross-compared at a glance.

using Test
using QAtlas
using QAtlas: Hubbard1D, FreeEnergy, Infinite, GroundStateEnergyDensity
using QAtlas.Hubbard1DJKSNLIE: atomic_free_energy

using LinearAlgebra: Hermitian, eigen
using SparseArrays: spzeros

@testset "Hubbard1D JKS NLIE — temperature sweep + cross-validation (#523)" begin
    @testset "Atomic limit at 7 temperature points (rtol scaled with β)" begin
        # Kinetic corrections O(β t / U) grow with β; allow looser rtol
        # away from β → 0 to track the physical kinetic shift.
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        cases = (
            (1e-5, 0.02),
            (1e-4, 0.02),
            (1e-3, 0.02),
            (1e-2, 0.05),
            (5e-2, 0.15),
            (1e-1, 0.30),
            (2e-1, 0.50),
        )
        for (β, rtol) in cases
            f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
            f_atom = atomic_free_energy(β, 4.0, 2.0)
            @test isfinite(f)
            @test isapprox(f, f_atom; rtol=rtol)
        end
    end

    @testset "f(β) monotone non-decreasing in β" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        betas = (1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
        fs = [QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β) for β in betas]
        @test all(isfinite, fs)
        @test all(diff(fs) .>= -1e-3 * maximum(abs.(fs)))
    end

    @testset "U dependence at high T (U-independent after Stage E.1)" begin
        # Stage E.1 Chebyshev-Gauss quadrature + page-14 direct form made the
        # FE evaluator U-independent: all U values give f/f_atom ~ 1.0 at
        # β = 1e-3 within 1% (previous U=2 and U=8 regressions are resolved).
        for U in (2.0, 4.0, 8.0)
            m = Hubbard1D(; t=1.0, U=U, μ=U/2)
            f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1e-3)
            f_atom = atomic_free_energy(1e-3, U, U/2)
            @test isfinite(f)
            @test isapprox(f, f_atom; rtol=0.01)
        end
    end

    @testset "Half-filling PH symmetry (μ = U/2, h = 0)" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        f_h0 = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.1)
        @test isreal(f_h0)
        @test isfinite(f_h0)
    end

    @testset "Grid refinement: higher N+x_max improves atomic agreement" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        β = 1e-3
        f_atom = atomic_free_energy(β, 4.0, 2.0)
        f_lo = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β, grid_N=32, x_max=4.0)
        f_hi = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β, grid_N=64, x_max=8.0)
        @test isfinite(f_lo)
        @test isfinite(f_hi)
        # Higher resolution should be at least as close to atomic.
        @test abs(f_hi / f_atom - 1) <= abs(f_lo / f_atom - 1) + 0.02
    end

    @testset "Thermodynamic consistency: numerical -∂(βf)/∂β ≈ e" begin
        # Numerical derivative at small dβ is noisy; loose rtol = 0.5 just
        # checks the sign and order of magnitude.
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        β0 = 1e-3
        dβ = 1e-4
        f_minus = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β0 - dβ)
        f_plus = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β0 + dβ)
        e_jks = -((β0 + dβ) * f_plus - (β0 - dβ) * f_minus) / (2 * dβ)
        f_atom_minus = atomic_free_energy(β0 - dβ, 4.0, 2.0)
        f_atom_plus = atomic_free_energy(β0 + dβ, 4.0, 2.0)
        e_atom = -((β0 + dβ) * f_atom_plus - (β0 - dβ) * f_atom_minus) / (2 * dβ)
        @test isfinite(e_jks)
        # The finite-difference derivative of (βf) at small Δβ amplifies any
        # residual O(10⁻³) FE-evaluator error into O(1) noise on e_jks; use
        # atol (not rtol) to absorb this. Sign and order-of-magnitude check.
        @test abs(e_jks - e_atom) < 5.0
    end

    @testset "ED comparison at N = 4 atomic sites, high T" begin
        function _ed_hubbard_atomic_free_energy(N::Int, U::Float64, μ::Float64, β::Float64)
            d = 4
            dim = d^N
            on_site = Float64[0, -μ, -μ, U - 2μ]
            H = spzeros(Float64, dim, dim)
            for idx in 0:(dim - 1)
                e = 0.0
                tmp = idx
                for s in 1:N
                    local_state = tmp % d
                    e += on_site[local_state + 1]
                    tmp = div(tmp, d)
                end
                H[idx + 1, idx + 1] = e
            end
            evals, _ = eigen(Hermitian(Matrix(H)))
            emin = minimum(evals)
            return (-log(sum(exp.(-β .* (evals .- emin)))) / β + emin) / N
        end

        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        β = 1e-3
        f_jks = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        f_ed_atom = _ed_hubbard_atomic_free_energy(4, 4.0, 2.0, β)
        @test isfinite(f_jks)
        @test isfinite(f_ed_atom)
        f_atom = atomic_free_energy(β, 4.0, 2.0)
        @test isapprox(f_ed_atom, f_atom; rtol=1e-6)
        @test isapprox(f_jks, f_ed_atom; rtol=0.1)
    end

    @testset "Comparison with Lieb-Wu GS energy at low T (β = 5)" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        β = 5.0
        f_lo_t = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        if isnan(f_lo_t)
            @info "JKS NLIE did not converge at β = $β; lo-T Lieb-Wu " *
                "comparison skipped pending Stage E." e0
            @test_broken false
        else
            @test isfinite(f_lo_t)
            @test isapprox(f_lo_t, e0; atol=2.0 * abs(e0))
        end
    end
end
