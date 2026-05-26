# ─────────────────────────────────────────────────────────────────────────────
# AKLT1D — structural / OBC-spectrum / identity guards (split, formerly testset 1 of 3)
#
# Split out of test/models/quantum/misc/test_aklt.jl (5.9 min on s02) so
# the three top-level testsets each run on their own shard. Helpers
# spin_ops, chain_hamiltonian, two_point, verify_profile_Ns come from
# test/util/{generic_ed,verify}.jl via runtests.jl ambient include.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

@testset "AKLT1D — structural / OBC-spectrum / identity guards" begin
    @testset "Constructor variants" begin
        @test AKLT1D().J == 1.0
        @test AKLT1D(; J=2.5).J == 2.5
        # J > 0 is required: every registered analytic observable assumes
        # the antiferromagnetic sign of the AKLT bond-projector form.
        @test_throws ArgumentError AKLT1D(; J=0.0)
        @test_throws ArgumentError AKLT1D(; J=-1.0)
        @test_throws ArgumentError AKLT1D(; J=(-eps()))
    end

    @testset "OBC ExactSpectrum shape (sorted, real, length 3^N)" begin
        m = AKLT1D(; J=1.0)
        for N in (2, 3, 4)
            λ = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
            @test length(λ) == 3^N
            @test issorted(λ)
            @test all(isreal, λ)
        end
    end

    @testset "4-fold OBC ground-state degeneracy (S_tot = 0 ⊕ 1)" begin
        # AKLT theorem: two free spin-1/2 edge modes → 4 = 1 + 3
        # degenerate ground states. Structural; not a single value.
        m = AKLT1D(; J=1.0)
        for N in (4, 6, 8)
            λ = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
            @test λ[4] - λ[1] < 1e-8                   # 4-fold manifold
            @test λ[5] - λ[4] > 0.05                    # well separated
        end
    end

    @testset "OBC ground-state energy scaling (analytic + bounds)" begin
        # E_0(N) = -(2/3)(N-1) J exact for OBC AKLT (every bond projector
        # gives zero). e_0/N → -2/3 with a 1/N edge correction.
        m = AKLT1D(; J=1.0)
        e_inf = -2 / 3
        for N in (4, 6, 8)
            E0 = QAtlas.fetch(m, ExactSpectrum(), OBC(N))[1]
            @test E0 ≈ -(2 / 3) * (N - 1) atol = 1e-10
            @test E0 / N > e_inf
            @test abs(E0 / N - e_inf) ≤ 2 / (3N)
        end
    end

    @testset "ZZCorrelation — symmetry, sign-alternation, ratio, J-independence" begin
        # Relational identities (not single-value cards):
        # r ↔ -r symmetry; alternating sign; ratio c_{r+1}/c_r = -1/3;
        # J-independent (VBS ground state J-invariant for J>0).
        m = AKLT1D(; J=1.0)
        for r in 1:5
            @test QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=(-r)) ≈
                QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=r) atol = 1e-14
        end
        c1 = QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=1)
        c2 = QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=2)
        c3 = QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=3)
        @test c1 < 0 && c2 > 0 && c3 < 0
        @test abs(c2 / c1) ≈ 1 / 3 atol = 1e-14
        @test abs(c3 / c2) ≈ 1 / 3 atol = 1e-14
        for J in (0.3, 1.0, 4.2)
            @test QAtlas.fetch(
                AKLT1D(; J=J), ZZCorrelation(; mode=:static), Infinite(); r=2
            ) ≈ QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=2) atol = 1e-14
        end
    end

    @testset "ZZStructureFactor — monotonicity, periodicity, J-independence" begin
        # Relational: 0 < S(q) < S(π) for 0 < q < π; 2π-periodic; even;
        # J-independent.
        m = AKLT1D(; J=1.0)
        Sπ = QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=π)
        for q in range(0.1, π - 0.1; length=12)
            S = QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=q)
            @test 0.0 < S < 2.0
            @test S < Sπ
        end
        @test QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=0.7) ≈
            QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=0.7 + 2π) atol = 1e-12
        @test QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=-0.7) ≈
            QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=0.7) atol = 1e-14
        @test QAtlas.fetch(AKLT1D(; J=2.6), ZZStructureFactor(), Infinite(); q=π) ≈ 2.0 atol =
            1e-14
    end

    @testset "Structure factor IS the Fourier transform of the correlation" begin
        # Cross-quantity identity: S_zz(q) = Σ_r e^{iqr} ⟨Sᶻ₀Sᶻ_r⟩.
        # Truncated lattice sum (geometric tail < 1e-12 by r = 40).
        m = AKLT1D(; J=1.0)
        for q in (0.3, 1.0, 2.0, π)
            Ssum = QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=0)
            for r in 1:40
                Ssum +=
                    2 *
                    cos(q * r) *
                    QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=r)
            end
            @test Ssum ≈ QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=q) atol = 1e-9
        end
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
