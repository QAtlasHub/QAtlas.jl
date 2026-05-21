# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: RMT (Wigner-Dyson β = 1, 2, 4) + Poisson universality.
#
# Migrated from pure-legacy @test to verify()-first (PR #449 phase B,
# zero-legacy end-state). Single-value pins (Atas mean ratio, Tracy-Widom
# F_β(0) Bornemann 2010 checkpoint, Poisson reference values) become
# verify() cards. Structural invariants (positivity, monotonicity,
# integral normalisations, asymptotic tails over loops) and
# error-throw paths stay raw @test — verify() is scalar-only.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, QuadGK

@testset "RMT: Wigner surmise — value at s = 0 and positivity" begin
    for β in (1, 2, 4)
        @test QAtlas.fetch(Universality(:RMT), WignerSurmise(); β=β, s=0.0) == 0.0
        for s in (0.1, 0.5, 1.0, 2.0, 5.0)
            v = QAtlas.fetch(Universality(:RMT), WignerSurmise(); β=β, s=s)
            @test v > 0
            @test isfinite(v)
        end
    end
end

@testset "RMT: Wigner surmise — small-s level repulsion P_β(s) ~ s^β" begin
    δ = 1e-4
    for β in (1, 2, 4)
        r =
            QAtlas.fetch(Universality(:RMT), WignerSurmise(); β=β, s=2δ) /
            QAtlas.fetch(Universality(:RMT), WignerSurmise(); β=β, s=δ)
        @test isapprox(r, 2.0^β; rtol=1e-3)
    end
end

@testset "RMT: Wigner surmise — normalisation ∫₀^∞ P_β(s) ds = 1" begin
    for β in (1, 2, 4)
        Pβ(s) = QAtlas.fetch(Universality(:RMT), WignerSurmise(); β=β, s=s)
        I, _ = quadgk(Pβ, 0.0, Inf; rtol=1e-10)
        @test isapprox(I, 1.0; atol=1e-8)
    end
end

@testset "RMT: Wigner surmise — first moment ∫₀^∞ s P_β(s) ds = 1" begin
    for β in (1, 2, 4)
        Pβ(s) = QAtlas.fetch(Universality(:RMT), WignerSurmise(); β=β, s=s)
        M, _ = quadgk(s -> s * Pβ(s), 0.0, Inf; rtol=1e-10)
        @test isapprox(M, 1.0; atol=1e-8)
    end
end

@testset "RMT: mean ratio ⟨r⟩ — Atas et al. 2013" begin
    for (β, r_lit) in ((1, 0.5307), (2, 0.5996), (4, 0.6744))
        verify(
            Universality(:RMT),
            MeanRatio(),
            Infinite();
            route=:literature_value,
            independent=r_lit,
            agree_within=1e-4,
            at=["β=$(β)"],
            refs=[
                "Atas-Bogomolny-Giraud-Roux 2013 (PRL 110, 084101): ⟨r⟩_β = $(r_lit) for Wigner-Dyson β=$(β)",
            ],
            fetch_kw=(; β=β),
        )
    end
end

@testset "RMT: Tracy-Widom F_β — Bornemann 2010 checkpoints at x = 0" begin
    for (β, F_lit) in ((1, 0.8319), (2, 0.9694), (4, 0.99966))
        verify(
            Universality(:RMT),
            TracyWidom(),
            Infinite();
            route=:literature_value,
            independent=F_lit,
            agree_within=1e-3,
            at=["β=$(β)", "x=0.0"],
            refs=[
                "Bornemann 2010 / Tracy-Widom 1994: F_$(β)(0) = $(F_lit) tabulated checkpoint",
            ],
            fetch_kw=(; β=β, x=0.0),
        )
    end
end

@testset "RMT: Tracy-Widom F_β — monotone non-decreasing on x ∈ [-5, 5]" begin
    grid = collect(range(-5.0, 5.0; step=0.25))
    for β in (1, 2, 4)
        vals = [QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=x) for x in grid]
        for i in 2:length(vals)
            @test vals[i] ≥ vals[i - 1] - 1e-12
        end
        @test all(0 ≤ v ≤ 1 for v in vals)
    end
end

@testset "RMT: Tracy-Widom F_β — tail behaviour" begin
    for β in (1, 2, 4)
        @test QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=-10.0) < 1e-3
        @test QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=-20.0) < 1e-30
        @test isapprox(
            QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=10.0), 1.0; atol=1e-3
        )
        @test isapprox(
            QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=50.0), 1.0; atol=1e-12
        )
    end
end

@testset "RMT: domain errors on unsupported β" begin
    for β_bad in (0, 3, 5, -1)
        @test_throws DomainError QAtlas.fetch(
            Universality(:RMT), WignerSurmise(); β=β_bad, s=1.0
        )
        @test_throws DomainError QAtlas.fetch(
            Universality(:RMT), TracyWidom(); β=β_bad, x=0.0
        )
        @test_throws DomainError QAtlas.fetch(Universality(:RMT), MeanRatio(); β=β_bad)
    end
end

@testset "Poisson: P(s) = exp(-s) — closed-form sweep (verify)" begin
    for s in (0.0, 0.5, 1.0, 2.0, 5.0)
        verify(
            Universality(:Poisson),
            WignerSurmise(),
            Infinite();
            route=:second_closed_form,
            independent=exp(-s),
            agree_within=1e-12,
            at=["s=$(s)"],
            refs=["Poisson level statistics: P(s) = exp(-s)"],
            fetch_kw=(; s=s),
        )
    end
    # NOTE: the standalone P(1) = 1/e card was redundant with the sweep
    # entry at s=1.0 (algebraically identical: exp(-1) = 1/e) and is
    # omitted to avoid a duplicate INVENTORY card.
end

@testset "Poisson: ⟨r⟩ = 2 log 2 - 1 (closed form, verify)" begin
    verify(
        Universality(:Poisson),
        MeanRatio(),
        Infinite();
        route=:second_closed_form,
        independent=2 * log(2) - 1,
        agree_within=1e-12,
        refs=["Poisson level-spacing mean ratio: ⟨r⟩ = 2 log 2 - 1 ≈ 0.386..."],
    )
end

@testset "Poisson: P(s) normalisation and first moment" begin
    P(s) = QAtlas.fetch(Universality(:Poisson), WignerSurmise(); s=s)
    I, _ = quadgk(P, 0.0, Inf; rtol=1e-10)
    M, _ = quadgk(s -> s * P(s), 0.0, Inf; rtol=1e-10)
    @test isapprox(I, 1.0; atol=1e-10)
    @test isapprox(M, 1.0; atol=1e-10)
end

@testset "RMT/Poisson: exported quantity types" begin
    @test WignerSurmise() isa QAtlas.AbstractQuantity
    @test TracyWidom() isa QAtlas.AbstractQuantity
    @test MeanRatio() isa QAtlas.AbstractQuantity
end

@testset "RMT: SpectralFormFactor — GUE late-time plateau K(τ→∞) = 1" begin
    verify(
        Universality(:RMT),
        SpectralFormFactor(),
        Infinite();
        route=:literature_value,
        independent=1.0,
        agree_within=0,
        refs=[
            "RMT GUE Heisenberg-time plateau: K(τ → ∞) = 1 (Mehta 1991; Bohigas-Giannoni-Schmit 1984)",
        ],
    )
    for τ in (2π, 2π + 1e-9, 10π, 100π, Inf)
        verify(
            Universality(:RMT),
            SpectralFormFactor(),
            Infinite();
            route=:literature_value,
            independent=1.0,
            agree_within=0,
            at=["τ=$(τ)"],
            refs=["GUE plateau holds for any τ ≥ 2π (Heisenberg time)"],
            fetch_kw=(; τ=τ),
        )
    end
    verify(
        Universality(:RMT),
        SpectralFormFactor(),
        Infinite();
        route=:literature_value,
        independent=1.0,
        agree_within=0,
        at=["ensemble=:GUE", "τ=Inf"],
        refs=["Explicit :GUE keyword identical to default"],
        fetch_kw=(; ensemble=:GUE, τ=Inf),
    )
end

@testset "RMT: SpectralFormFactor — Phase 2 deferrals raise DomainError" begin
    for τ_bad in (0.0, 0.5, 1.0, π, 2π - 1e-6)
        @test_throws DomainError QAtlas.fetch(
            Universality(:RMT), SpectralFormFactor(); τ=τ_bad
        )
    end
    @test_throws DomainError QAtlas.fetch(
        Universality(:RMT), SpectralFormFactor(); ensemble=:GOE
    )
    @test_throws DomainError QAtlas.fetch(
        Universality(:RMT), SpectralFormFactor(); ensemble=:GSE
    )
    @test_throws DomainError QAtlas.fetch(
        Universality(:RMT), SpectralFormFactor(); ensemble=:CUE
    )
end

@testset "RMT: SpectralFormFactor — exported quantity type" begin
    @test SpectralFormFactor() isa QAtlas.AbstractQuantity
end
