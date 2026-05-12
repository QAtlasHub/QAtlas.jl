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
