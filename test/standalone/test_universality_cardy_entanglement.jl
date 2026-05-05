# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Calabrese–Cardy entanglement entropy at the
# Universality{C} level (issue #149).
#
# Targeted run:
#   julia --project=test test/standalone/test_universality_cardy_entanglement.jl
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality :: Cardy entanglement entropy" begin

    # ── Ising c = 1/2: closed-form value at ℓ = 4, L = 8 ─────────────────────
    @testset "Ising d=2 (c=1/2) PBC closed form" begin
        c = 1 / 2
        ℓ, L = 4.0, 8.0
        expected = (c / 3) * log((L / π) * sin(π * ℓ / L))
        got = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=ℓ, L=L)
        @test got ≈ expected atol = 1e-10
    end

    # ── PBC vs OBC coefficient ratio ─────────────────────────────────────────
    #
    # The bare prefactor ratio is exactly 2:1 (c/3 vs c/6).  The OBC log
    # argument has an extra factor of 2 inside (image doubling); after
    # subtracting the (c/6) log 2 shift the OBC value is exactly half of
    # the PBC value.
    @testset "PBC vs OBC: OBC = (PBC/2) up to log-2 image shift" begin
        ℓ, L = 50.0, 100.0
        c = 1 / 2
        Spbc = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=ℓ, L=L)
        Sobc = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), OBC(); ℓ=ℓ, L=L)
        @test (Sobc - (c / 6) * log(2)) ≈ Spbc / 2 atol = 1e-12
    end

    # ── Infinite limit ───────────────────────────────────────────────────────
    @testset "Infinite (PBC scaling) c=1/2" begin
        c = 1 / 2
        ℓ = 10.0
        expected = (c / 3) * log(ℓ)
        got = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), Infinite(); ℓ=ℓ)
        @test got ≈ expected atol = 1e-12
    end

    # ── Maximum at ℓ = L/2 ; UV divergences at ℓ=0, ℓ=L ──────────────────────
    @testset "PBC maximum at ℓ=L/2 and divergences at endpoints" begin
        L = 64.0
        Smid = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=L / 2, L=L)
        Soff1 = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=L / 4, L=L)
        Soff2 = QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=3 * L / 8, L=L
        )
        @test Smid ≥ Soff1
        @test Smid ≥ Soff2
        # Endpoints: ℓ = 0 and ℓ = L give log(0) → -Inf
        Send_left = QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=0.0, L=L
        )
        Send_right = QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=L, L=L
        )
        @test Send_left == -Inf
        @test Send_right == -Inf
    end

    # ── Rényi α=2 reduces to Cardy with c → c · 3/4 ──────────────────────────
    #
    # Substitution c -> c · (1 + 1/α) / 2 at α = 2 gives c_eff = (3/4) c.
    # The PBC Rényi prefactor is c_eff/3 = (3/4)·(c/3), i.e. 3/4 of the
    # von Neumann result with the same log argument.
    @testset "Rényi α=2 reduces to (3/4)·VN" begin
        ℓ, L = 4.0, 8.0
        S_vn = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=ℓ, L=L)
        S_r2 = QAtlas.fetch(Universality(:Ising), RenyiEntropy(2.0), PBC(); ℓ=ℓ, L=L)
        @test S_r2 ≈ (3 // 4) * S_vn atol = 1e-12

        S_vn_obc = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), OBC(); ℓ=ℓ, L=L)
        S_r2_obc = QAtlas.fetch(Universality(:Ising), RenyiEntropy(2.0), OBC(); ℓ=ℓ, L=L)
        @test S_r2_obc ≈ (3 // 4) * S_vn_obc atol = 1e-12

        S_vn_inf = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), Infinite(); ℓ=ℓ)
        S_r2_inf = QAtlas.fetch(Universality(:Ising), RenyiEntropy(2.0), Infinite(); ℓ=ℓ)
        @test S_r2_inf ≈ (3 // 4) * S_vn_inf atol = 1e-12
    end

    # ── Other classes: Potts3 (c=4/5), Potts4 (c=1), XY (c=1), Heisenberg (c=1)
    @testset "Potts3 / Potts4 / XY / Heisenberg central charges" begin
        @test QAtlas.fetch(Universality(:Potts3), CentralCharge()) == 4 // 5
        @test QAtlas.fetch(Universality(:Potts4), CentralCharge()) == 1 // 1
        @test QAtlas.fetch(Universality(:XY), CentralCharge()) == 1 // 1
        @test QAtlas.fetch(Universality(:Ising), CentralCharge()) == 1 // 2
        @test QAtlas.fetch(Universality(:Heisenberg), CentralCharge(); d=1) == 1 // 1

        ℓ, L = 4.0, 8.0
        S_potts3 = QAtlas.fetch(Universality(:Potts3), VonNeumannEntropy(), PBC(); ℓ=ℓ, L=L)
        @test S_potts3 ≈ ((4 / 5) / 3) * log((L / π) * sin(π * ℓ / L)) atol = 1e-12

        S_potts4_inf = QAtlas.fetch(
            Universality(:Potts4), VonNeumannEntropy(), Infinite(); ℓ=10.0
        )
        @test S_potts4_inf ≈ (1 / 3) * log(10.0) atol = 1e-12

        S_xy_obc = QAtlas.fetch(Universality(:XY), VonNeumannEntropy(), OBC(); ℓ=ℓ, L=L)
        @test S_xy_obc ≈ (1 / 6) * log((2 * L / π) * sin(π * ℓ / L)) atol = 1e-12
    end

    # ── Class without CentralCharge defined → ErrorException ────────────────
    @testset "Classes without CentralCharge raise ErrorException" begin
        @test_throws ErrorException QAtlas.fetch(
            Universality(:KPZ), VonNeumannEntropy(), PBC(); ℓ=4.0, L=8.0
        )
        @test_throws ErrorException QAtlas.fetch(
            Universality(:Percolation), VonNeumannEntropy(), Infinite(); ℓ=10.0
        )
    end

    # ── d ≥ 3 Ising / XY / Heisenberg d≥2: not a 1+1D CFT, error out ────────
    @testset "non-1+1D classes have no central charge" begin
        @test_throws ErrorException QAtlas.fetch(Universality(:Ising), CentralCharge(); d=3)
        @test_throws ErrorException QAtlas.fetch(Universality(:XY), CentralCharge(); d=3)
        @test_throws ErrorException QAtlas.fetch(Universality(:Heisenberg), CentralCharge(); d=2)
        @test_throws ErrorException QAtlas.fetch(Universality(:Heisenberg), CentralCharge(); d=3)
    end

    # ── Argument validation ─────────────────────────────────────────────────
    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=-1.0, L=8.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=10.0, L=8.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=4.0, L=-1.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), Infinite(); ℓ=0.0
        )
    end
end
