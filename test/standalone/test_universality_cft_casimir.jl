# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: CFT Casimir / finite-size ground-state energy correction
# (Cardy 1986, Blöte–Cardy–Nightingale 1986, Affleck 1986)
#
# Phase 1 of issue #150 — only the universal 1/L correction is exercised.
# The `ConformalTower` quantity is Phase 2 / out of scope here.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality: CFT Casimir correction (Cardy 1986)" begin
    # ─── Ising c = 1/2 PBC, exact closed-form check ─────────────────────────
    @testset "Ising c=1/2 PBC reference value" begin
        # E_0^PBC 1/L term = -π c v / (6 L)
        # at c = 1/2, v = 2.0, L = 16.0 → -π · (1/2) · 2 / (6 · 16) = -π/96
        val = QAtlas.fetch(
            Universality(:Ising), CasimirEnergyCorrection(), PBC(); L=16.0, v=2.0
        )
        @test val ≈ -π / 96 atol = 1e-12
    end

    # ─── Ising c = 1/2 OBC ──────────────────────────────────────────────────
    @testset "Ising c=1/2 OBC reference value" begin
        val = QAtlas.fetch(
            Universality(:Ising), CasimirEnergyCorrection(), OBC(); L=16.0, v=2.0
        )
        # at c = 1/2, v = 2.0, L = 16.0 → -π · (1/2) · 2 / (24 · 16) = -π/384
        @test val ≈ -π / 384 atol = 1e-12
    end

    # ─── PBC : OBC ratio = 4, class-independent (kinematic CFT result) ──────
    @testset "PBC : OBC ratio = 4, class-independent" begin
        for class in (:Ising, :Potts3, :Potts4, :XY, :Heisenberg)
            pbc = QAtlas.fetch(
                Universality(class), CasimirEnergyCorrection(), PBC(); L=20.0, v=1.7
            )
            obc = QAtlas.fetch(
                Universality(class), CasimirEnergyCorrection(), OBC(); L=20.0, v=1.7
            )
            @test pbc ≈ 4 * obc atol = 1e-12
        end
    end

    # ─── L → ∞ limit → 0 (1/L scaling) ──────────────────────────────────────
    @testset "Thermodynamic limit: 1/L scaling, value -> 0" begin
        v = 2.0
        for L in (1e3, 1e4, 1e5, 1e6)
            pbc = QAtlas.fetch(
                Universality(:Ising), CasimirEnergyCorrection(), PBC(); L=L, v=v
            )
            # 1/L scaling: pbc * L is L-independent
            @test pbc * L ≈ -π * (1 // 2) * v / 6 atol = 1e-12
        end
        # vanishes at very large L
        huge = QAtlas.fetch(
            Universality(:Ising), CasimirEnergyCorrection(), PBC(); L=1e10, v=2.0
        )
        @test abs(huge) < 1e-9
    end

    # ─── Multi-class numeric values (Potts3 c=4/5, Potts4 c=1, XY/Heis c=1) ─
    @testset "Per-class central charge values" begin
        # Potts3 c = 4/5
        p3 = QAtlas.fetch(
            Universality(:Potts3), CasimirEnergyCorrection(), PBC(); L=10.0, v=1.0
        )
        @test p3 ≈ -π * (4 // 5) / 60 atol = 1e-12  # -π · (4/5) · 1 / (6·10)
        # Potts4 c = 1
        p4 = QAtlas.fetch(
            Universality(:Potts4), CasimirEnergyCorrection(), PBC(); L=10.0, v=1.0
        )
        @test p4 ≈ -π / 60 atol = 1e-12
        # XY c = 1
        xy = QAtlas.fetch(
            Universality(:XY), CasimirEnergyCorrection(), PBC(); L=10.0, v=1.0
        )
        @test xy ≈ -π / 60 atol = 1e-12
        # Heisenberg c = 1
        hb = QAtlas.fetch(
            Universality(:Heisenberg), CasimirEnergyCorrection(), PBC(); L=10.0, v=1.0
        )
        @test hb ≈ -π / 60 atol = 1e-12
    end

    # ─── Unsupported classes raise ErrorException ───────────────────────────
    @testset "Unsupported classes raise ErrorException" begin
        # KPZ: non-equilibrium, no CFT central charge
        @test_throws ErrorException QAtlas.fetch(
            Universality(:KPZ), CasimirEnergyCorrection(), PBC(); L=16.0, v=2.0
        )
        # Percolation: non-unitary logarithmic CFT, c = 0 not via Cardy
        @test_throws ErrorException QAtlas.fetch(
            Universality(:Percolation), CasimirEnergyCorrection(), OBC(); L=16.0, v=2.0
        )
        # Unknown class
        @test_throws ErrorException QAtlas.fetch(
            Universality(:Bogus), CasimirEnergyCorrection(), PBC(); L=16.0, v=2.0
        )
    end

    # ─── Domain checks: L, v must be positive ───────────────────────────────
    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), CasimirEnergyCorrection(), PBC(); L=-1.0, v=2.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), CasimirEnergyCorrection(), PBC(); L=16.0, v=0.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), CasimirEnergyCorrection(), OBC(); L=0.0, v=2.0
        )
    end
end
