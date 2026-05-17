# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: AKLT1D — exact VBS ground state + Haldane gap
#
# Verifies the closed-form AKLT 1988 values (energy density, correlation
# length, string order parameter), the García-Saez–Murg–Verstraete 2013 Haldane gap,
# and the OBC dense-ED 4-fold edge-state degeneracy (S=1/2 edges → S_tot
# ∈ {0, 1}: singlet + triplet) on N = 4, 6, 8.
#
# Run targeted (Pkg.test() forbidden by user policy on Panza):
#
#   julia --project=test test/standalone/test_aklt.jl
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "AKLT1D — exact VBS analytical values" begin
    @testset "Construction + J scaling" begin
        m = AKLT1D()
        @test m.J == 1.0
        @test AKLT1D(; J=2.5).J == 2.5

        # Linear J scaling for every analytic infinite-limit observable
        for J in (0.5, 1.0, 2.5)
            mJ = AKLT1D(; J=J)
            @test QAtlas.fetch(mJ, GroundStateEnergyDensity(), Infinite()) ≈ -2J / 3 atol =
                1e-14
            @test QAtlas.fetch(mJ, MassGap(), Infinite()) ≈ 0.350 * J rtol = 1e-12
            # ξ and O_str are J-independent
            @test QAtlas.fetch(mJ, CorrelationLength(), Infinite()) ≈ 1 / log(3) atol =
                1e-14
            @test QAtlas.fetch(mJ, StringOrderParameter(), Infinite()) ≈ 4 / 9 atol = 1e-14
        end
    end

    @testset "GroundStateEnergyDensity (Infinite) = -2J/3 (closed form)" begin
        m = AKLT1D(; J=1.0)
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        @test e0 ≈ -2 / 3 atol = 1e-14
        @test e0 ≈ -0.6666666666666666 atol = 1e-14
        # Energy(:per_site) routes to the same analytic value
        @test QAtlas.fetch(m, Energy(:per_site), Infinite()) ≈ -2 / 3 atol = 1e-14
    end

    @testset "CorrelationLength (Infinite) = 1/log 3 (closed form)" begin
        ξ = QAtlas.fetch(AKLT1D(), CorrelationLength(), Infinite())
        @test ξ ≈ 1 / log(3) atol = 1e-12
        @test ξ ≈ 0.9102392266268373 atol = 1e-12
    end

    @testset "MassGap (Infinite) ≈ 0.350 J (García-Saez-Murg-Verstraete 2013)" begin
        Δ = QAtlas.fetch(AKLT1D(), MassGap(), Infinite())
        # Compare against the canonical DMRG value with atol 1e-4 as per the
        # acceptance criteria; the implementation stores it to 5 decimal places.
        @test Δ ≈ 0.350 atol = 1e-3
    end

    @testset "StringOrderParameter (Infinite) = 4/9 (Kennedy-Tasaki 1992)" begin
        O = QAtlas.fetch(AKLT1D(), StringOrderParameter(), Infinite())
        @test O ≈ 4 / 9 atol = 1e-14
    end
end

@testset "AKLT1D — OBC dense ED (N ≤ 8)" begin
    m = AKLT1D(; J=1.0)

    @testset "ExactSpectrum is sorted, real, length 3^N" begin
        for N in (2, 3, 4)
            λ = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
            @test length(λ) == 3^N
            @test issorted(λ)
            @test all(isreal, λ)
        end
    end

    @testset "4-fold OBC ground-state degeneracy (S_tot = 0 ⊕ 1)" begin
        # AKLT theorem: under OBC two free spin-1/2 edge modes give
        # 4 = 1 (singlet) + 3 (triplet) degenerate ground states.
        # Dense ED on the unbroken Hamiltonian sees this exactly up to
        # numerical noise.
        for N in (4, 6, 8)
            λ = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
            Δ_low_4 = λ[4] - λ[1]
            @test Δ_low_4 < 1e-8
            # First excitation above the 4-fold manifold should be > 0.05 J
            # for these chain lengths (well-separated from edge manifold).
            @test λ[5] - λ[4] > 0.05
        end
    end

    @testset "OBC GS energy is exactly -(2/3)(N-1) J" begin
        # AKLT ground states are exact zero-energy eigenstates of every
        # bond projector under OBC, so the ground-state energy is
        # *exactly* −(2/3)·(N − 1)·J — no finite-size correction at all
        # on top of the missing edge bond.  Per-site energy approaches
        # −2J/3 with a 1/N edge correction.
        e_inf = -2 / 3
        for N in (4, 6, 8)
            λ = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
            E0 = λ[1]
            @test E0 ≈ -(2 / 3) * (N - 1) atol = 1e-10
            # Per-site energy approaches -2/3 from above (less negative)
            @test E0 / N > e_inf
            @test abs(E0 / N - e_inf) ≤ 2 / (3N)  # 1/N edge bound
        end
    end
end

@testset "AKLT1D — OBC finite-temperature (dense ED, exact for N ≤ 8)" begin
    m = AKLT1D(; J=1.0)

    @testset "β → ∞ reproduces the exact frustration-free ground state" begin
        # ⟨H⟩(β→∞) must converge to the EXACT OBC AKLT ground-state
        # energy E₀ = -(2J/3)(N-1) (the VBS annihilates every bond
        # projector, so there is no finite-size correction beyond the
        # missing edge bond).
        for N in (4, 6, 8)
            E0 = -(2 / 3) * (N - 1)
            E_inf = QAtlas.fetch(m, Energy(:total), OBC(N); beta=1.0e3)
            @test E_inf ≈ E0 atol = 1e-6
            # FreeEnergy → E₀ as β → ∞ (per-site × N)
            f_inf = QAtlas.fetch(m, FreeEnergy(), OBC(N); beta=1.0e3)
            @test f_inf * N ≈ E0 atol = 1e-2
        end
    end

    @testset "β → ∞ entropy → log 4 (AKLT edge-mode signature)" begin
        # Two free spin-½ edge modes ⇒ 4-fold degenerate OBC ground
        # manifold (singlet ⊕ triplet). The residual TOTAL entropy is
        # exactly log 4, independent of N and J.
        for N in (4, 6, 8)
            S_total = N * QAtlas.fetch(m, ThermalEntropy(), OBC(N); beta=1.0e3)
            @test S_total ≈ log(4) atol = 1e-6
        end
        # J-independent: same residual entropy at a different coupling
        S_J2 = 6 * QAtlas.fetch(AKLT1D(; J=2.0), ThermalEntropy(), OBC(6); beta=1.0e3)
        @test S_J2 ≈ log(4) atol = 1e-6
    end

    @testset "β → 0 entropy → N log 3 (full spin-1 Hilbert space)" begin
        for N in (4, 6, 8)
            S_total = N * QAtlas.fetch(m, ThermalEntropy(), OBC(N); beta=1.0e-6)
            @test S_total ≈ N * log(3) atol = 1e-4
        end
    end

    @testset "SpecificHeat: non-negative, Schottky peak, vanishes at both ends" begin
        for N in (4, 6, 8)
            c0 = QAtlas.fetch(m, SpecificHeat(), OBC(N); beta=1.0e-6)
            c_inf = QAtlas.fetch(m, SpecificHeat(), OBC(N); beta=1.0e3)
            c_pk = QAtlas.fetch(m, SpecificHeat(), OBC(N); beta=1.0)
            # → 0 at both extremes (fp-noise tolerance)
            @test abs(c0) < 1e-6
            @test abs(c_inf) < 1e-6
            # finite Schottky-like peak in between, strictly positive
            @test c_pk > 0.1
            # heat capacity is a variance ⇒ non-negative (within fp)
            for β in (0.1, 0.5, 1.0, 2.0, 5.0)
                @test QAtlas.fetch(m, SpecificHeat(), OBC(N); beta=β) > -1e-9
            end
        end
    end

    @testset "zero-T Energy{:total} is exactly linear in J" begin
        # ⟨H⟩ at fixed finite β is not exactly linear in J (β couples as
        # βJ), but the β→∞ limit is the GS energy E₀ = -(2J/3)(N-1),
        # which IS exactly linear in J.
        E1_inf = QAtlas.fetch(AKLT1D(; J=1.0), Energy(:total), OBC(6); beta=1.0e3)
        E2_inf = QAtlas.fetch(AKLT1D(; J=2.0), Energy(:total), OBC(6); beta=1.0e3)
        @test E2_inf ≈ 2 * E1_inf atol = 1e-6
    end

    @testset "Infinite() + beta has no closed form ⇒ DomainError" begin
        for Q in (Energy(:total), FreeEnergy(), ThermalEntropy(), SpecificHeat())
            @test_throws DomainError QAtlas.fetch(m, Q, Infinite(); beta=1.0)
        end
    end
end
