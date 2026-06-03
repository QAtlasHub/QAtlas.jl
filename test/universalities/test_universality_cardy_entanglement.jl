# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Calabrese-Cardy entanglement entropy at the
# Universality{C} level (issue #149).
#
# Migrated from pure-legacy @test to verify()-first (PR #449 phase B,
# zero-legacy end-state). Value pins (closed forms + literature CentralCharge +
# Renyi → VN delegation) become verify() cards; structural inequalities and
# -Inf edge pins stay raw because verify() requires finite scalar best.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality :: Cardy entanglement entropy" begin
    # ── Ising c = 1/2 PBC: closed-form value at ℓ = 4, L = 8 ────────────────
    let c = 1 / 2, ℓ = 4.0, L = 8.0
        verify(
            Universality(:Ising),
            VonNeumannEntropy(),
            PBC();
            route=:second_closed_form,
            independent=(c / 3) * log((L / π) * sin(π * ℓ / L)),
            agree_within=1e-10,
            refs=[
                "Calabrese-Cardy 2004: S_PBC(ℓ; L) = (c/3) log[(L/π) sin(πℓ/L)] at c=1/2"
            ],
            fetch_kw=(; ℓ=ℓ, L=L),
        )
    end

    # ── PBC vs OBC log-2 image shift ────────────────────────────────────────
    let c = 1 / 2, ℓ = 50.0, L = 100.0
        Spbc = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=ℓ, L=L)
        verify(
            Universality(:Ising),
            VonNeumannEntropy(),
            OBC();
            route=:delegation_invariant,
            independent=Spbc / 2 + (c / 6) * log(2),
            agree_within=1e-12,
            refs=["Calabrese-Cardy 2004 image construction: S_OBC = S_PBC/2 + (c/6) log 2"],
            fetch_kw=(; ℓ=ℓ, L=L),
        )
    end

    # ── Infinite (PBC scaling) c=1/2 ────────────────────────────────────────
    let c = 1 / 2, ℓ = 10.0
        verify(
            Universality(:Ising),
            VonNeumannEntropy(),
            Infinite();
            route=:second_closed_form,
            independent=(c / 3) * log(ℓ),
            agree_within=1e-12,
            refs=["Calabrese-Cardy 2004: S_∞(ℓ) = (c/3) log ℓ + const at c=1/2"],
            fetch_kw=(; ℓ=ℓ),
        )
    end

    # ── Structural: PBC max at ℓ=L/2 + UV-divergent endpoints
    # (NOT single-value hub pins; kept raw @test) ────────────────────────────
    @testset "PBC maximum at ℓ=L/2 and divergences at endpoints (structural)" begin
        L = 64.0
        Smid = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=L / 2, L=L)
        Soff1 = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=L / 4, L=L)
        Soff2 = QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=3 * L / 8, L=L
        )
        @test Smid ≥ Soff1
        @test Smid ≥ Soff2
        Send_left = QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=0.0, L=L
        )
        Send_right = QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=L, L=L
        )
        @test Send_left == -Inf
        @test Send_right == -Inf
    end

    # ── Rényi α=2 reduces to (3/4)·VN — delegation_invariant ────────────────
    let ℓ = 4.0, L = 8.0
        S_vn_pbc = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), PBC(); ℓ=ℓ, L=L)
        verify(
            Universality(:Ising),
            RenyiEntropy(2.0),
            PBC();
            route=:delegation_invariant,
            independent=(3 // 4) * S_vn_pbc,
            agree_within=1e-12,
            refs=[
                "Calabrese-Cardy 2004 Rényi substitution c → c·(1 + 1/α)/2 ⇒ S_R(α=2) = (3/4) S_VN with same log argument",
            ],
            fetch_kw=(; ℓ=ℓ, L=L),
        )

        S_vn_obc = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), OBC(); ℓ=ℓ, L=L)
        verify(
            Universality(:Ising),
            RenyiEntropy(2.0),
            OBC();
            route=:delegation_invariant,
            independent=(3 // 4) * S_vn_obc,
            agree_within=1e-12,
            refs=["Calabrese-Cardy Rényi α=2: S_R^OBC = (3/4) S_VN^OBC"],
            fetch_kw=(; ℓ=ℓ, L=L),
        )

        S_vn_inf = QAtlas.fetch(Universality(:Ising), VonNeumannEntropy(), Infinite(); ℓ=ℓ)
        verify(
            Universality(:Ising),
            RenyiEntropy(2.0),
            Infinite();
            route=:delegation_invariant,
            independent=(3 // 4) * S_vn_inf,
            agree_within=1e-12,
            refs=["Calabrese-Cardy Rényi α=2 (Infinite): S_R = (3/4) S_VN"],
            fetch_kw=(; ℓ=ℓ),
        )
    end

    # ── Other-class CentralCharge values + non-Ising SF closed forms ────────
    @testset "Potts3 / Potts4 / XY / Heisenberg central charges + SF" begin
        # CentralCharge literature pins kept as raw @test here: the verify()
        # 3-arg dispatch fetch(Universality{X}, CentralCharge, ::Infinite)
        # is added by PR #451 (on main) but not yet on `next` (this PR's
        # base). Verify-card coverage of these constants lives in PR #451's
        # standalone test_universality_central_charge_lit.jl. Once #451
        # lands on next, this block can be migrated to verify() too.
        @test QAtlas.fetch(Universality(:Potts3), CentralCharge()) == 4 // 5
        @test QAtlas.fetch(Universality(:Potts4), CentralCharge()) == 1 // 1
        @test QAtlas.fetch(Universality(:XY), CentralCharge()) == 1 // 1
        @test QAtlas.fetch(Universality(:Ising), CentralCharge()) == 1 // 2
        @test QAtlas.fetch(Universality(:Heisenberg), CentralCharge(); d=1) == 1 // 1

        # Closed-form SF cross-checks on the Cardy formula.
        let ℓ = 4.0, L = 8.0
            verify(
                Universality(:Potts3),
                VonNeumannEntropy(),
                PBC();
                route=:second_closed_form,
                independent=((4 / 5) / 3) * log((L / π) * sin(π * ℓ / L)),
                agree_within=1e-12,
                refs=["Calabrese-Cardy at c=4/5 (Potts3) PBC"],
                fetch_kw=(; ℓ=ℓ, L=L),
            )
            verify(
                Universality(:Potts4),
                VonNeumannEntropy(),
                Infinite();
                route=:second_closed_form,
                independent=(1 / 3) * log(10.0),
                agree_within=1e-12,
                refs=["Calabrese-Cardy at c=1 (Potts4) Infinite"],
                fetch_kw=(; ℓ=10.0),
            )
            verify(
                Universality(:XY),
                VonNeumannEntropy(),
                OBC();
                route=:second_closed_form,
                independent=(1 / 6) * log((2 * L / π) * sin(π * ℓ / L)),
                agree_within=1e-12,
                refs=[
                    "Calabrese-Cardy OBC at c=1 (XY): factor 1/6, image-doubled log argument",
                ],
                fetch_kw=(; ℓ=ℓ, L=L),
            )
        end
    end

    # ── Error-path guards (raw @test_throws — verify() doesn't model error
    # outcomes) ──────────────────────────────────────────────────────────────
    @testset "Classes without CentralCharge raise ErrorException" begin
        @test_throws ErrorException QAtlas.fetch(
            Universality(:KPZ), VonNeumannEntropy(), PBC(); ℓ=4.0, L=8.0
        )
        @test_throws ErrorException QAtlas.fetch(
            Universality(:Percolation), VonNeumannEntropy(), Infinite(); ℓ=10.0
        )
    end

    @testset "non-1+1D classes have no central charge" begin
        @test_throws ErrorException QAtlas.fetch(Universality(:Ising), CentralCharge(); d=3)
        @test_throws ErrorException QAtlas.fetch(Universality(:XY), CentralCharge(); d=3)
        @test_throws ErrorException QAtlas.fetch(
            Universality(:Heisenberg), CentralCharge(); d=2
        )
        @test_throws ErrorException QAtlas.fetch(
            Universality(:Heisenberg), CentralCharge(); d=3
        )
    end

    @testset "Infinite finite-T VN (CC sinh form, #580)" begin
        # New beta kwarg: S(ℓ, β) = (c/3) log[(β/π) sinh(πℓ/β)] across all
        # 1+1D-CFT universality classes registered. β → ∞ recovers T=0.
        # Each class has a canonical dimension that yields a 1+1D CFT
        # at the universality level; Heisenberg uses d=1 (1D chain).
        d_for_class = Dict(:Ising=>2, :XY=>2, :Heisenberg=>1, :Potts3=>2, :Potts4=>2)
        for class in (:Ising, :XY, :Heisenberg, :Potts3, :Potts4)
            d = d_for_class[class]
            c = float(QAtlas.fetch(Universality(class), CentralCharge(); d=d))
            for ℓ in (5.0, 25.0), β in (1.0, 10.0, 100.0)
                S = QAtlas.fetch(
                    Universality(class), VonNeumannEntropy(), Infinite(); ℓ=ℓ, beta=β
                )
                expected = (c / 3) * log((β / π) * sinh(π * ℓ / β))
                @test S ≈ expected atol = 1e-12
            end
            # β → ∞ recovery
            S_T0 = QAtlas.fetch(
                Universality(class), VonNeumannEntropy(), Infinite(); ℓ=10.0, beta=Inf
            )
            S_lowT = QAtlas.fetch(
                Universality(class), VonNeumannEntropy(), Infinite(); ℓ=10.0, beta=1.0e6
            )
            @test S_lowT ≈ S_T0 rtol = 1e-3
        end
    end

    @testset "Infinite finite-T Renyi (CC sinh form with c → c·(1+1/α)/2, #580)" begin
        d_for_class = Dict(:Ising=>2, :Heisenberg=>1)
        for class in (:Ising, :Heisenberg)
            c = float(
                QAtlas.fetch(Universality(class), CentralCharge(); d=d_for_class[class])
            )
            for α in (0.5, 2.0, 5.0), ℓ in (10.0, 30.0), β in (2.0, 20.0)
                S = QAtlas.fetch(
                    Universality(class), RenyiEntropy(α), Infinite(); ℓ=ℓ, beta=β
                )
                c_eff = c * (1 + 1 / α) / 2
                expected = (c_eff / 3) * log((β / π) * sinh(π * ℓ / β))
                @test S ≈ expected atol = 1e-12
            end
        end
    end

    @testset "Infinite finite-T argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), VonNeumannEntropy(), Infinite(); ℓ=10.0, beta=-1.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), RenyiEntropy(2.0), Infinite(); ℓ=10.0, beta=0.0
        )
    end

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
