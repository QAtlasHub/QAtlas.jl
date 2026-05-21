# ─────────────────────────────────────────────────────────────────────────────
# Test: Majumdar–Ghosh chain — structural / error / identity / relational
# guards that verify() architecturally cannot express (@test_throws,
# constructor invariants, registry sanity, deprecation-alias behaviour,
# SU(2) cross-checks, strict bounds).
#
# Split out of test_majumdar_ghosh.jl (10.4 min on s01) — this part is
# trivially fast; the heavy ED sweeps live in test_majumdar_ghosh_verify_ed_*.jl,
# the literature / closed-form pins in test_majumdar_ghosh_verify_lit.jl.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using Logging: with_logger, NullLogger

@testset "MajumdarGhosh — structural / error / identity guards" begin
    @testset "Constructor variants agree" begin
        @test MajumdarGhosh(1.0).J == MajumdarGhosh(; J=1.0).J
    end

    @testset "GSED PBC — odd N rejected, N-kwarg form callable" begin
        m = MajumdarGhosh(; J=1.0)
        @test_throws DomainError QAtlas.fetch(m, GroundStateEnergyDensity(), PBC(5))
        # N-kwarg API form reachable (value covered by verify cards below).
        @test QAtlas.fetch(m, GroundStateEnergyDensity(), PBC(); N=8) isa Real
    end

    @testset "MassGap method dispatch & deprecation" begin
        m = MajumdarGhosh(; J=1.0)
        Δ_num = QAtlas.fetch(m, MassGap(), Infinite(); method=:numerical)
        Δ_trimer = QAtlas.fetch(m, MassGap(), Infinite(); method=:trimer_bound)
        # Default routes to :numerical (no fetch_kw == same path).
        @test QAtlas.fetch(m, MassGap(), Infinite()) == Δ_num
        # Strict relational: SS-1981 trimer-sector bound exceeds the
        # actual gap (cannot be expressed by a single-value verify card).
        @test Δ_trimer > Δ_num
        # Legacy :lower_bound alias resolves to J/4 with a one-shot @warn
        # (deprecation behaviour; not a verify-card concern).
        Δ_legacy = with_logger(NullLogger()) do
            QAtlas.fetch(m, MassGap(), Infinite(); method=:lower_bound)
        end
        @test Δ_legacy ≈ 0.25 atol = 1e-14
        # Unsupported method symbols raise DomainError.
        @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); method=:dmrg_strict)
        @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); method=:bogus)
    end

    @testset "Registry knows about MajumdarGhosh" begin
        rows = QAtlas.implementation_status(MajumdarGhosh)
        @test !isempty(rows)
        quantities = unique(r.quantity for r in rows)
        @test GroundStateEnergyDensity in quantities
        @test MassGap in quantities
    end

    @testset "SpinGap — strict bounds, DomainError, SU(2) identity" begin
        # Strict bounds (relational; not single-value verify cards):
        #   Δ < J/4  (Shastry-Sutherland 1981 trimer-sector bound)
        #   Δ > 0.117 J  (Magnus 1991 strict absolute-gap bound)
        Δ = QAtlas.fetch(MajumdarGhosh(), SpinGap(), Infinite())
        @test Δ < 0.25
        @test Δ > 0.117

        # SpinGap rejects J ≤ 0.
        @test_throws DomainError QAtlas.fetch(MajumdarGhosh(), SpinGap(), Infinite(); J=0.0)
        @test_throws DomainError QAtlas.fetch(
            MajumdarGhosh(), SpinGap(), Infinite(); J=-1.5
        )

        # SU(2) identity: MG is SU(2)-symmetric, so the spectral gap
        # equals the spin gap (S=0 → S=1 triplet excitation).
        for J in (0.5, 1.0, 3.0)
            m = MajumdarGhosh(; J=J)
            @test QAtlas.fetch(m, SpinGap(), Infinite()) ≈
                QAtlas.fetch(m, MassGap(), Infinite(); method=:numerical)
        end
    end
end
