# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: RMT (Wigner-Dyson β = 1, 2, 4) + Poisson universality.
#
# Phase 1 of issue #151:
#   * Wigner surmise (closed-form) for β = 1, 2, 4
#   * Tracy-Widom F_β via Bornemann 2010 table + tail asymptotics
#   * Atas-Bogomolny-Giraud-Roux ⟨r⟩ for β = 1, 2, 4 and Poisson
#   * Poisson level statistics (P(s) = exp(-s))
# Painlevé-II direct integrator deferred to Phase 2.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, QuadGK

# ─── Wigner surmise ──────────────────────────────────────────────────────────

@testset "RMT: Wigner surmise — value at s = 0 and positivity" begin
    for β in (1, 2, 4)
        # P_β(0) = 0 by construction (level repulsion, exact).
        @test QAtlas.fetch(Universality(:RMT), WignerSurmise(); β=β, s=0.0) == 0.0
        # Positive on (0, ∞).
        for s in (0.1, 0.5, 1.0, 2.0, 5.0)
            v = QAtlas.fetch(Universality(:RMT), WignerSurmise(); β=β, s=s)
            @test v > 0
            @test isfinite(v)
        end
    end
end

@testset "RMT: Wigner surmise — small-s level repulsion P_β(s) ~ s^β" begin
    # Compare ratios: P_β(2δ) / P_β(δ) should approach 2^β as δ → 0.
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

# ─── Mean ratio ⟨r⟩ ──────────────────────────────────────────────────────────

@testset "RMT: mean ratio ⟨r⟩ — Atas et al. 2013" begin
    @test isapprox(QAtlas.fetch(Universality(:RMT), MeanRatio(); β=1), 0.5307; atol=1e-4)
    @test isapprox(QAtlas.fetch(Universality(:RMT), MeanRatio(); β=2), 0.5996; atol=1e-4)
    @test isapprox(QAtlas.fetch(Universality(:RMT), MeanRatio(); β=4), 0.6744; atol=1e-4)
end

# ─── Tracy-Widom F_β ─────────────────────────────────────────────────────────

@testset "RMT: Tracy-Widom F_β — Bornemann 2010 checkpoints at x = 0" begin
    @test isapprox(
        QAtlas.fetch(Universality(:RMT), TracyWidom(); β=1, x=0.0), 0.8319; atol=1e-3
    )
    @test isapprox(
        QAtlas.fetch(Universality(:RMT), TracyWidom(); β=2, x=0.0), 0.9694; atol=1e-3
    )
    @test isapprox(
        QAtlas.fetch(Universality(:RMT), TracyWidom(); β=4, x=0.0), 0.99966; atol=1e-3
    )
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
        # Left tail → 0 for x → -∞.
        @test QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=-10.0) < 1e-3
        @test QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=-20.0) < 1e-30
        # Right tail → 1 for x → +∞.
        @test isapprox(
            QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=10.0), 1.0; atol=1e-3
        )
        @test isapprox(
            QAtlas.fetch(Universality(:RMT), TracyWidom(); β=β, x=50.0), 1.0; atol=1e-12
        )
    end
end

# ─── DomainError on unsupported β ────────────────────────────────────────────

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

# ─── Poisson universality ────────────────────────────────────────────────────

@testset "Poisson: P(s) = exp(-s)" begin
    for s in (0.0, 0.5, 1.0, 2.0, 5.0)
        @test isapprox(
            QAtlas.fetch(Universality(:Poisson), WignerSurmise(); s=s), exp(-s); atol=1e-12
        )
    end
    # Pinned reference: P(1) = 1/e.
    @test isapprox(
        QAtlas.fetch(Universality(:Poisson), WignerSurmise(); s=1.0), 1 / exp(1); atol=1e-12
    )
end

@testset "Poisson: ⟨r⟩ = 2 log 2 - 1" begin
    @test isapprox(
        QAtlas.fetch(Universality(:Poisson), MeanRatio()), 2 * log(2) - 1; atol=1e-12
    )
end

@testset "Poisson: P(s) normalisation and first moment" begin
    P(s) = QAtlas.fetch(Universality(:Poisson), WignerSurmise(); s=s)
    I, _ = quadgk(P, 0.0, Inf; rtol=1e-10)
    M, _ = quadgk(s -> s * P(s), 0.0, Inf; rtol=1e-10)
    @test isapprox(I, 1.0; atol=1e-10)
    @test isapprox(M, 1.0; atol=1e-10)
end

# ─── Type-export sanity ──────────────────────────────────────────────────────

@testset "RMT/Poisson: exported quantity types" begin
    @test WignerSurmise() isa QAtlas.AbstractQuantity
    @test TracyWidom() isa QAtlas.AbstractQuantity
    @test MeanRatio() isa QAtlas.AbstractQuantity
end
