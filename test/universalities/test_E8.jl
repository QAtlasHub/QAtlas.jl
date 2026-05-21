# ─────────────────────────────────────────────────────────────────────────────
# E8 spectrum (Zamolodchikov 1989) — value pins via verify(), structural
# guards + legacy-alias deprecation shims via raw @test.
#
# Migrated from pure-legacy @test to verify()-first (PR #449 phase B,
# zero-legacy end-state). E8Spectrum returns a Vector{Float64} of 8 mass
# ratios; verify() takes a scalar so subject_extract = m -> m[i] projects
# each component. Coldea et al. (Science 2010) measured the m₂/m₁ = φ
# ratio in CoNb₂O₆ — that experimental cross-check anchors the literature
# routes here.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "E8 Spectrum Logic" begin
    # ── Structural sanity (length, ordering, positivity) — raw @test
    # (multi-element invariants of the Vector, not single-value hub pins)
    @testset "Structure" begin
        masses = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        @test length(masses) == 8
        @test masses[1] == 1.0
        @test all(masses .> 0)
        @test issorted(masses)
    end

    # ── Golden ratio m₂/m₁ = φ — Coldea et al. Science 327, 177 (2010) ──────
    verify(
        E8(),
        E8Spectrum(),
        Infinite();
        route=:literature_value,
        independent=(1 + sqrt(5)) / 2,
        agree_within=1e-15,
        refs=[
            "Coldea et al., Science 327, 177 (2010): m₂/m₁ = φ measured in CoNb₂O₆ Ising chain at h=h_c",
        ],
        subject_extract=m -> m[2],
    )

    # ── Closed-form mass spectrum (Zamolodchikov 1989, Delfino 2004) ────────
    @testset "Closed-form expressions (Zamolodchikov 1989, Delfino 2004)" begin
        ϕ = 2 * cos(π / 5)
        closedforms = [
            ("m1=1", 1.0, m -> m[1]),
            ("m2=ϕ", ϕ, m -> m[2]),
            ("m3", 2 * cos(π / 30), m -> m[3]),
            ("m4=2c(7π/30)ϕ", 2 * cos(7π / 30) * ϕ, m -> m[4]),
            ("m5=2c(2π/15)ϕ", 2 * cos(2π / 15) * ϕ, m -> m[5]),
            ("m6=2c(π/30)ϕ", 2 * cos(π / 30) * ϕ, m -> m[6]),
            ("m7=2c(7π/30)ϕ²", 2 * cos(7π / 30) * ϕ^2, m -> m[7]),
            ("m8=2c(2π/15)ϕ²", 2 * cos(2π / 15) * ϕ^2, m -> m[8]),
        ]
        for (label, expected, extract) in closedforms
            verify(
                E8(),
                E8Spectrum(),
                Infinite();
                route=:second_closed_form,
                independent=expected,
                agree_within=1e-14,
                at=[label],
                refs=[
                    "Zamolodchikov 1989 + Delfino 2004: Perron-Frobenius eigenvector of the E₈ Cartan matrix gives the 8-mass spectrum in closed form; see docs/src/calc/e8-mass-spectrum-derivation.md",
                ],
                subject_extract=extract,
            )
        end
    end

    # ── Tabulated decimal values (Delfino 2004, eq. 4.14) ───────────────────
    @testset "Tabulated decimal values (Delfino 2004, eq. 4.14)" begin
        expected = [
            1.000000, 1.618034, 1.989044, 2.404867, 2.956295, 3.218340, 3.891157, 4.783386
        ]
        for i in 1:8
            verify(
                E8(),
                E8Spectrum(),
                Infinite();
                route=:literature_value,
                independent=expected[i],
                agree_within=1e-5,
                at=["m$(i)"],
                refs=[
                    "Delfino 2004 eq. (4.14) / Zamolodchikov 1989 Table 2: six-digit tabulated value of m_$(i)",
                ],
                subject_extract=m -> m[i],
            )
        end
    end

    # ── Fusion-rule multiplicative identities (Zamolodchikov 1989 §5) ───────
    @testset "Fusion-rule multiplicative identities" begin
        m_ref = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        for (i, j, k) in ((2, 3, 6), (2, 4, 7), (2, 5, 8))
            verify(
                E8(),
                E8Spectrum(),
                Infinite();
                route=:delegation_invariant,
                independent=m_ref[i] * m_ref[j],
                agree_within=1e-14,
                at=["fusion $(i)×$(j)→$(k)"],
                refs=[
                    "Zamolodchikov 1989 §5 / Delfino 2004 Table 2: on-shell bootstrap fusion a×b→c forces m_c = m_a m_b for these triples",
                ],
                subject_extract=m -> m[k],
            )
        end
    end

    # ── Legacy Symbol-dispatch shim — kept raw (tests deprecation @info path)
    @testset "Aliases" begin
        expected = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        @test @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(:E8, :mass_ratios) ==
            expected
        @test @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(:E8, :E8_masses) ==
            expected
        @test @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(:E8, :mass_ratio) ==
            expected
    end

    @testset "Type Stability" begin
        @test @inferred(QAtlas.fetch(Model(:E8), Quantity(:E8_spectrum), Infinite())) isa
            Vector{Float64}
    end
end
