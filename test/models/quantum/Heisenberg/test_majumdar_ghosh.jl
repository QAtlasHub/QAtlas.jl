# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Majumdar–Ghosh chain (S = 1/2 J₁–J₂ at J₂/J₁ = 1/2)
#
# Verify the exact dimer-product ground-state energy density E₀/N = −3J/8
# (Majumdar–Ghosh 1969), the analytical Shastry–Sutherland (1981) gap
# lower bound Δ ≥ J/4, and the White–Affleck (1996) DMRG numerical gap
# Δ ≈ 0.234 J.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using Logging: with_logger, NullLogger

@testset "MajumdarGhosh: ground state and gap" begin
    @testset "GroundStateEnergyDensity (Infinite) — exact dimer GS" begin
        # J = 1: closed-form -3/8 (size-independent dimer-product state)
        m = MajumdarGhosh()
        e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        @test e0 == -3 / 8

        # J scaling: E₀/N is linear in J
        for Jval in (0.5, 2.0, 3.7)
            mJ = MajumdarGhosh(; J=Jval)
            @test QAtlas.fetch(mJ, GroundStateEnergyDensity(), Infinite()) == -3 * Jval / 8
        end

        # Constructor variants must agree
        @test MajumdarGhosh(1.0).J == MajumdarGhosh(; J=1.0).J
    end

    @testset "GroundStateEnergyDensity (PBC) — size-independent" begin
        m = MajumdarGhosh(; J=1.0)
        for N in (4, 8, 12, 16)
            e0 = QAtlas.fetch(m, GroundStateEnergyDensity(), PBC(N))
            @test e0 == -3 / 8
        end

        # Odd N is rejected (no consistent dimer covering on a ring).
        @test_throws DomainError QAtlas.fetch(m, GroundStateEnergyDensity(), PBC(5))

        # N kwarg form also works.
        @test QAtlas.fetch(m, GroundStateEnergyDensity(), PBC(); N=8) == -3 / 8
    end

    @testset "MassGap (Infinite) — :numerical default + :trimer_bound + legacy alias" begin
        m = MajumdarGhosh(; J=1.0)

        # Default method is now :numerical (the actual DMRG ground-state-
        # to-first-excited gap, ≈ 0.234 J).  Switched from the previous
        # :lower_bound default after PR review noted that the SS 1981
        # J/4 value exceeds the true gap and is best read as a sector-
        # specific bound.
        Δ_default = QAtlas.fetch(m, MassGap(), Infinite())
        Δ_num_explicit = QAtlas.fetch(m, MassGap(), Infinite(); method=:numerical)
        @test Δ_default == Δ_num_explicit
        @test Δ_default ≈ 0.234 atol = 1e-14

        # Shastry-Sutherland trimer-sector bound.
        Δ_trimer = QAtlas.fetch(m, MassGap(), Infinite(); method=:trimer_bound)
        @test Δ_trimer ≈ 0.25 atol = 1e-14
        @test Δ_trimer > Δ_num_explicit  # SS bound numerically exceeds the actual gap

        # Legacy :lower_bound alias still resolves to J/4 (with a one-shot
        # @warn pointing callers to :trimer_bound).  Wrapped to silence
        # the deprecation warning during regression testing.
        Δ_legacy = with_logger(NullLogger()) do
            QAtlas.fetch(m, MassGap(), Infinite(); method=:lower_bound)
        end
        @test Δ_legacy ≈ 0.25 atol = 1e-14

        # J scaling for both physical methods.
        for Jval in (0.5, 2.0, 3.7)
            mJ = MajumdarGhosh(; J=Jval)
            @test QAtlas.fetch(mJ, MassGap(), Infinite(); method=:trimer_bound) == Jval / 4
            @test QAtlas.fetch(mJ, MassGap(), Infinite(); method=:numerical) ≈ 0.234 * Jval atol =
                1e-14
        end

        # Unsupported method symbols raise DomainError.
        @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); method=:dmrg_strict)
        @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); method=:bogus)
    end

    @testset "Registry rows" begin
        # Basic sanity that the registry knows about MajumdarGhosh.
        rows = QAtlas.implementation_status(MajumdarGhosh)
        @test !isempty(rows)
        quantities = unique(r.quantity for r in rows)
        @test GroundStateEnergyDensity in quantities
        @test MassGap in quantities
    end

    @testset "MajumdarGhosh — SpinGap Δ ≈ 0.234 J (Phase 2, White-Affleck 1996)" begin
        Δ = QAtlas.fetch(MajumdarGhosh(), SpinGap(), Infinite())
        @test Δ ≈ 0.234
        # Linear scaling with J
        @test QAtlas.fetch(MajumdarGhosh(; J=3.0), SpinGap(), Infinite()) ≈ 3 * 0.234
        @test QAtlas.fetch(MajumdarGhosh(; J=0.5), SpinGap(), Infinite()) ≈ 0.5 * 0.234
        # SpinGap < Shastry-Sutherland trimer bound J/4 (DMRG < analytical sector bound, as expected)
        @test Δ < 0.25
        # SpinGap > strict Magnus 1991 absolute-gap bound 0.117 J
        @test Δ > 0.117
    end

    @testset "MajumdarGhosh — SpinGap rejects J ≤ 0 (Phase 2)" begin
        @test_throws DomainError QAtlas.fetch(MajumdarGhosh(), SpinGap(), Infinite(); J=0.0)
        @test_throws DomainError QAtlas.fetch(
            MajumdarGhosh(), SpinGap(), Infinite(); J=-1.5
        )
    end

    @testset "MajumdarGhosh — third-pass: SpinGap == MassGap(method=:numerical) (SU(2) cross-check)" begin
        # MG is SU(2)-symmetric; the lowest excitation is a triplet, so the
        # spectral gap (MassGap, :numerical) equals the spin gap (S=0 → S=1).
        for J in (0.5, 1.0, 3.0)
            m = MajumdarGhosh(; J=J)
            @test QAtlas.fetch(m, SpinGap(), Infinite()) ≈
                QAtlas.fetch(m, MassGap(), Infinite(); method=:numerical)
        end
    end
end
