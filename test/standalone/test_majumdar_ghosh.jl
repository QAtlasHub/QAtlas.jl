# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Majumdar–Ghosh chain (S = 1/2 J₁–J₂ at J₂/J₁ = 1/2)
#
# Verify the exact dimer-product ground-state energy density E₀/N = −3J/8
# (Majumdar–Ghosh 1969), the analytical Shastry–Sutherland (1981) gap
# lower bound Δ ≥ J/4, and the White–Affleck (1996) DMRG numerical gap
# Δ ≈ 0.234 J.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

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

    @testset "MassGap (Infinite) — :lower_bound vs :numerical" begin
        m = MajumdarGhosh(; J=1.0)

        # Default method is the analytical Shastry–Sutherland bound.
        Δ_lb_default = QAtlas.fetch(m, MassGap(), Infinite())
        Δ_lb_explicit = QAtlas.fetch(m, MassGap(), Infinite(); method=:lower_bound)
        @test Δ_lb_default == Δ_lb_explicit
        @test Δ_lb_default ≈ 0.25 atol = 1e-14

        # White–Affleck DMRG numerical-exact gap.
        Δ_num = QAtlas.fetch(m, MassGap(), Infinite(); method=:numerical)
        @test Δ_num ≈ 0.234 atol = 0.01

        # J scaling for both methods.
        for Jval in (0.5, 2.0, 3.7)
            mJ = MajumdarGhosh(; J=Jval)
            @test QAtlas.fetch(mJ, MassGap(), Infinite(); method=:lower_bound) == Jval / 4
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
end
