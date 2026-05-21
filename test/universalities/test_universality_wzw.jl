# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: WZW SU(2)_k universality classes
#
# Migrated from pure-legacy @test to verify()-first (PR #449 phase B,
# zero-legacy end-state). Sugawara CentralCharge c(k) = 3k/(k+2) and
# primary conformal weights h_j = j(j+1)/(k+2) become verify() cards
# (literature_value / second_closed_form); construction errors and
# Rational{Int} type-stability remain raw @test (verify() can not model
# error-throwing dispatch or type-trait checks).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

# ═══════════════════ Construction & validation (raw) ═════════════════════════

@testset "WZWSU2: construction validation" begin
    @test_throws DomainError WZWSU2(0)
    @test_throws DomainError WZWSU2(-1)
    @test WZWSU2(1).k == 1
    @test WZWSU2(7).k == 7
end

# ═══════════════════ Central charge — Sugawara closed form ═══════════════════

@testset "WZWSU2: CentralCharge (Sugawara)" begin
    for (k, c_lit) in ((1, 1 // 1), (2, 3 // 2), (3, 9 // 5), (4, 2 // 1), (10, 5 // 2))
        verify(
            WZWSU2(k),
            CentralCharge(),
            Infinite();
            route=:second_closed_form,
            independent=c_lit,
            agree_within=0,
            at=["k=$(k)"],
            refs=[
                "Knizhnik-Zamolodchikov 1984 (Sugawara): c(k) = 3k/(k+2) for SU(2)_k WZW"
            ],
        )
    end
    # Type stability (multi-element trait check, not single-value pin → raw).
    for k in (1, 2, 3, 7, 25)
        @test QAtlas.fetch(WZWSU2(k), CentralCharge()) isa Rational{Int}
    end
end

# ═══════════════════ Primary conformal weights — closed form ═════════════════

@testset "WZWSU2: k=1 primary weights (j ∈ {0, 1/2})" begin
    verify(
        WZWSU2(1),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=0 // 1,
        agree_within=0,
        refs=["Knizhnik-Zamolodchikov 1984: h_j = j(j+1)/(k+2); j=0, k=1 ⇒ 0"],
        fetch_kw=(; j=0),
    )
    verify(
        WZWSU2(1),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=1 // 4,
        agree_within=0,
        refs=["j=1/2, k=1 ⇒ (1/2)(3/2)/3 = 1/4 (free-boson at SU(2) radius)"],
        fetch_kw=(; j=1 // 2),
    )
    @test_throws DomainError QAtlas.fetch(WZWSU2(1), ConformalWeights(); j=1)
end

@testset "WZWSU2: k=2 primary weights (j ∈ {0, 1/2, 1})" begin
    verify(
        WZWSU2(2),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=0 // 1,
        agree_within=0,
        refs=["j=0, k=2 ⇒ 0"],
        fetch_kw=(; j=0),
    )
    verify(
        WZWSU2(2),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=3 // 16,
        agree_within=0,
        refs=["j=1/2, k=2 ⇒ (1/2)(3/2)/4 = 3/16"],
        fetch_kw=(; j=1 // 2),
    )
    verify(
        WZWSU2(2),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=1 // 2,
        agree_within=0,
        refs=["j=1, k=2 ⇒ 1·2/4 = 1/2"],
        fetch_kw=(; j=1),
    )
    @test_throws DomainError QAtlas.fetch(WZWSU2(2), ConformalWeights(); j=3 // 2)
end

@testset "WZWSU2: k=3 primary weights" begin
    verify(
        WZWSU2(3),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=0 // 1,
        agree_within=0,
        refs=["j=0, k=3 ⇒ 0"],
        fetch_kw=(; j=0),
    )
    verify(
        WZWSU2(3),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=3 // 20,
        agree_within=0,
        refs=["j=1/2, k=3 ⇒ (1/2)(3/2)/5 = 3/20"],
        fetch_kw=(; j=1 // 2),
    )
    verify(
        WZWSU2(3),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=2 // 5,
        agree_within=0,
        refs=["j=1, k=3 ⇒ 1·2/5 = 2/5"],
        fetch_kw=(; j=1),
    )
    verify(
        WZWSU2(3),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=3 // 4,
        agree_within=0,
        refs=["j=3/2, k=3 ⇒ (3/2)(5/2)/5 = 15/20 = 3/4"],
        fetch_kw=(; j=3 // 2),
    )
    @test_throws DomainError QAtlas.fetch(WZWSU2(3), ConformalWeights(); j=2)
end

@testset "WZWSU2: k=4 primary weights" begin
    verify(
        WZWSU2(4),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=1 // 1,
        agree_within=0,
        refs=["j=2, k=4 ⇒ 2·3/6 = 1"],
        fetch_kw=(; j=2),
    )
    verify(
        WZWSU2(4),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=1 // 3,
        agree_within=0,
        refs=["j=1, k=4 ⇒ 1·2/6 = 1/3"],
        fetch_kw=(; j=1),
    )
    @test_throws DomainError QAtlas.fetch(WZWSU2(4), ConformalWeights(); j=5 // 2)
end

# ═══════════════════ Validation: half-integer / non-negative (raw) ═══════════

@testset "WZWSU2: invalid j → DomainError" begin
    w = WZWSU2(4)
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=-1 // 2)
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=-1)
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=1 // 3)
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=2 // 3)
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=0.5)
end

# ═══════════════════ Type stability (raw) ═══════════════════════════════════

@testset "WZWSU2: ConformalWeights returns Rational{Int}" begin
    @test QAtlas.fetch(WZWSU2(1), ConformalWeights(); j=0) isa Rational{Int}
    @test QAtlas.fetch(WZWSU2(1), ConformalWeights(); j=1 // 2) isa Rational{Int}
    @test QAtlas.fetch(WZWSU2(7), ConformalWeights(); j=3) isa Rational{Int}
end
