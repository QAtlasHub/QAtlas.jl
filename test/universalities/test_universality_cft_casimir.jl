# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: CFT Casimir / finite-size ground-state energy correction
# (Cardy 1986, Blöte–Cardy–Nightingale 1986, Affleck 1986)
#
# Migrated from pure-legacy @test to verify()-first (PR #449 phase B,
# zero-legacy end-state). Value pins now run through verify() so each
# (class, bc, L, v) tuple becomes a structural INVENTORY card.
# Error-path @test_throws stay raw — verify() cannot represent
# dispatches that intentionally throw.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality: CFT Casimir correction (Cardy 1986)" begin
    # ─── Ising c = 1/2 PBC, exact closed form ───────────────────────────────
    verify(
        Universality(:Ising),
        CasimirEnergyCorrection(),
        PBC();
        route=:second_closed_form,
        independent=-π / 96,
        agree_within=1e-12,
        refs=[
            "Cardy 1986: E_0^PBC universal 1/L correction = -π c v / (6 L); Ising c=1/2, v=2.0, L=16 ⇒ -π/96",
        ],
        fetch_kw=(; L=16.0, v=2.0),
    )

    # ─── Ising c = 1/2 OBC ──────────────────────────────────────────────────
    verify(
        Universality(:Ising),
        CasimirEnergyCorrection(),
        OBC();
        route=:second_closed_form,
        independent=-π / 384,
        agree_within=1e-12,
        refs=[
            "Cardy 1986 / Blöte-Cardy-Nightingale 1986: E_0^OBC = -π c v / (24 L); Ising c=1/2, v=2.0, L=16 ⇒ -π/384",
        ],
        fetch_kw=(; L=16.0, v=2.0),
    )

    # ─── PBC : OBC ratio = 4, class-independent (Cardy 1986 kinematic) ──────
    @testset "PBC : OBC ratio = 4, class-independent" begin
        for class in (:Ising, :Potts3, :Potts4, :XY, :Heisenberg)
            verify(
                Universality(class),
                CasimirEnergyCorrection(),
                PBC();
                route=:delegation_invariant,
                independent=4 * QAtlas.fetch(
                    Universality(class), CasimirEnergyCorrection(), OBC(); L=20.0, v=1.7
                ),
                agree_within=1e-12,
                at=["class=$(class)"],
                refs=[
                    "Cardy 1986 kinematic CFT result: E_0^PBC / E_0^OBC = (1/6)/(1/24) = 4, independent of central charge / boundary conditions / class",
                ],
                fetch_kw=(; L=20.0, v=1.7),
            )
        end
    end

    # ─── Thermodynamic limit: 1/L scaling ───────────────────────────────────
    @testset "Thermodynamic limit: 1/L scaling, value -> 0" begin
        v = 2.0
        for L in (1e3, 1e4, 1e5, 1e6)
            verify(
                Universality(:Ising),
                CasimirEnergyCorrection(),
                PBC();
                route=:second_closed_form,
                independent=-π * (1 // 2) * v / 6 / L,    # = -π c v / (6 L) at c=1/2
                agree_within=1e-12,
                at=["L=$(L)"],
                refs=[
                    "Cardy 1986 1/L scaling: pbc·L = -π·c·v/6 is L-independent ⇒ pbc = -π·c·v/(6L)",
                ],
                fetch_kw=(; L=L, v=v),
            )
        end
        # vanishes at very large L (no verify card — just a tail bound)
        huge = QAtlas.fetch(
            Universality(:Ising), CasimirEnergyCorrection(), PBC(); L=1e10, v=2.0
        )
        @test abs(huge) < 1e-9
    end

    # ─── Multi-class numeric values (Potts3 c=4/5, Potts4 c=1, XY/Heis c=1) ─
    @testset "Per-class central charge values" begin
        # Potts3 c = 4/5: -π · (4/5) · 1 / (6 · 10) = -π · 4 / 300 = -π · (4//5) / 60
        verify(
            Universality(:Potts3),
            CasimirEnergyCorrection(),
            PBC();
            route=:literature_value,
            independent=-π * (4 // 5) / 60,
            agree_within=1e-12,
            refs=[
                "Cardy 1986 + Dotsenko 1984: 2D 3-state Potts c=4/5 ⇒ E_0^PBC/L = -π·c·v/(6L) = -π·(4/5)/60 at v=1, L=10",
            ],
            fetch_kw=(; L=10.0, v=1.0),
        )
        # Potts4 c = 1
        verify(
            Universality(:Potts4),
            CasimirEnergyCorrection(),
            PBC();
            route=:literature_value,
            independent=-π / 60,
            agree_within=1e-12,
            refs=[
                "Cardy 1986 + DFMS §12.3: 2D 4-state Potts (marginal compact boson) c=1 ⇒ E_0^PBC/L = -π/60 at v=1, L=10",
            ],
            fetch_kw=(; L=10.0, v=1.0),
        )
        # XY c = 1
        verify(
            Universality(:XY),
            CasimirEnergyCorrection(),
            PBC();
            route=:literature_value,
            independent=-π / 60,
            agree_within=1e-12,
            refs=[
                "Cardy 1986 + Kosterlitz 1974: 2D XY (BKT free compact boson) c=1 ⇒ E_0^PBC/L = -π/60 at v=1, L=10",
            ],
            fetch_kw=(; L=10.0, v=1.0),
        )
        # Heisenberg c = 1
        verify(
            Universality(:Heisenberg),
            CasimirEnergyCorrection(),
            PBC();
            route=:literature_value,
            independent=-π / 60,
            agree_within=1e-12,
            refs=[
                "Cardy 1986 + Affleck-Haldane 1987: spin-1/2 Heisenberg chain (SU(2)_1 WZW) c=1 ⇒ E_0^PBC/L = -π/60 at v=1, L=10",
            ],
            fetch_kw=(; L=10.0, v=1.0),
        )
    end

    # ─── Error-path guards (kept raw @test_throws — verify() doesn't model
    # error-throwing dispatch outcomes) ────────────────────────────────────
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
