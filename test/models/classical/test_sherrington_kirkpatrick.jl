# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: SherringtonKirkpatrick — T_c = J.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "SherringtonKirkpatrick — T_c = J" begin
    for J in (0.5, 1.0, 2.5)
        @test QAtlas.fetch(
            SherringtonKirkpatrick(; J=J), CriticalTemperature(), Infinite()
        ) ≈ J
    end
end

@testset "SherringtonKirkpatrick — J ≤ 0 returns 0" begin
    @test QAtlas.fetch(
        SherringtonKirkpatrick(; J=0.0), CriticalTemperature(), Infinite()
    ) == 0.0
    @test QAtlas.fetch(
        SherringtonKirkpatrick(; J=-1.5), CriticalTemperature(), Infinite()
    ) == 0.0
end

@testset "SherringtonKirkpatrick — Parisi T=0 ground-state energy density (Phase 2)" begin
    # Default J=1 → e_0 ≈ -0.7631667
    e0 = QAtlas.fetch(SherringtonKirkpatrick(), Energy{:per_site}(), Infinite())
    # Crisanti-Rizzo 2002 quote precision ±1e-5 — match literature, do NOT use
    # default √eps tolerance which would break on legit future refinements.
    @test isapprox(e0, -0.7631667; atol=1e-5)
    @test e0 < 0  # ground-state is negative for spin-glass mean-field
    # Rigorous lower bound: the SK Hamiltonian is bounded below by the
    # annealed / RS bound; full-RSB lifts this further but e_0 > -1 is
    # a robust sanity guard against accidental sign flip / overflow.
    @test e0 > -1.0
    # Scales linearly with J
    e0_3 = QAtlas.fetch(SherringtonKirkpatrick(; J=3.0), Energy{:per_site}(), Infinite())
    @test e0_3 ≈ 3 * e0
    # Identifies the Parisi/full-RSB value within Crisanti-Rizzo error bar
    @test isapprox(e0, -0.7631667; atol=1e-5)
end

@testset "SherringtonKirkpatrick — Energy rejects J ≤ 0 (Phase 2)" begin
    m = SherringtonKirkpatrick(; J=1.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); J=0.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); J=-1.0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "SherringtonKirkpatrick — verification cards" begin
    # Spin-glass transition Tc = J (mean-field saddle point)
    for J in (1.0, 2.0)
        verify(
            SherringtonKirkpatrick(; J=J),
            CriticalTemperature(),
            Infinite();
            route=:second_closed_form,
            independent=Float64(J),
            agree_within=1e-12,
            refs=["SK mean-field: spin-glass transition at Tc = J"],
        )
    end

    # T=0 Parisi full-RSB ground-state energy density (literature numeric)
    verify(
        SherringtonKirkpatrick(; J=1.0),
        Energy(:per_site),
        Infinite();
        route=:literature_value,
        independent=-0.7631667,
        agree_within=1e-4,
        refs=["Crisanti-Rizzo 2002 / Parisi full-RSB: e0 ≈ -0.7631667 (J=1)"],
    )
end
# ── additional verification cards (#381 batch 3) ─────────────────────────
@testset "SherringtonKirkpatrick — CriticalTemperature (#381 batch 3)" begin
    # Sherrington-Kirkpatrick spin glass mean-field critical temperature
    # T_c = J (in J-units, k_B = 1) — replica-symmetric and full-RSB
    # transition (Sherrington-Kirkpatrick 1975; Parisi 1979).
    for J in (0.5, 1.0, 2.0)
        verify(
            SherringtonKirkpatrick(; J=J),
            CriticalTemperature(),
            Infinite();
            route=:second_closed_form,
            independent=Float64(J),
            agree_within=1e-12,
            refs=[
                "Sherrington-Kirkpatrick 1975: SK spin-glass T_c = J in mean field (k_B = 1)",
            ],
        )
    end
end
