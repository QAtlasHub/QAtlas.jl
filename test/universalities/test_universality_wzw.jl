# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: WZW SU(2)_k universality classes
#
# Verifies the Sugawara central charge c(k) = 3k / (k + 2) and the
# primary conformal weights h_j = j(j+1) / (k+2) as exact
# `Rational{Int}` values, plus the validation rules for the spin
# parameter j (half-integer, 0 ≤ j ≤ k/2).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

# ═══════════════════ Construction & validation ════════════════════════════════

@testset "WZWSU2: construction validation" begin
    @test_throws DomainError WZWSU2(0)
    @test_throws DomainError WZWSU2(-1)
    @test WZWSU2(1).k == 1
    @test WZWSU2(7).k == 7
end

# ═══════════════════ Central charge — exact rationals ═════════════════════════

@testset "WZWSU2: CentralCharge (Sugawara)" begin
    @test QAtlas.fetch(WZWSU2(1), CentralCharge()) == 1 // 1     # free boson at SU(2) radius
    @test QAtlas.fetch(WZWSU2(2), CentralCharge()) == 3 // 2     # Ising × Majorana
    @test QAtlas.fetch(WZWSU2(3), CentralCharge()) == 9 // 5
    @test QAtlas.fetch(WZWSU2(4), CentralCharge()) == 2 // 1
    @test QAtlas.fetch(WZWSU2(10), CentralCharge()) == 30 // 12  # 5 // 2
    # Type stability.
    for k in (1, 2, 3, 7, 25)
        @test QAtlas.fetch(WZWSU2(k), CentralCharge()) isa Rational{Int}
    end
end

# ═══════════════════ k = 1: Heisenberg / free-boson cross-check ═══════════════

@testset "WZWSU2: k=1 primary weights (j ∈ {0, 1/2})" begin
    w = WZWSU2(1)
    @test QAtlas.fetch(w, ConformalWeights(); j=0) == 0 // 1
    @test QAtlas.fetch(w, ConformalWeights(); j=1 // 2) == 1 // 4
    # Above max spin should error.
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=1)
end

# ═══════════════════ k = 2: full primary list ═════════════════════════════════

@testset "WZWSU2: k=2 primary weights (j ∈ {0, 1/2, 1})" begin
    w = WZWSU2(2)
    @test QAtlas.fetch(w, ConformalWeights(); j=0) == 0 // 1
    @test QAtlas.fetch(w, ConformalWeights(); j=1 // 2) == 3 // 16
    @test QAtlas.fetch(w, ConformalWeights(); j=1) == 1 // 2
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=3 // 2)
end

# ═══════════════════ Higher k spot-checks ═════════════════════════════════════

@testset "WZWSU2: k=3 primary weights" begin
    w = WZWSU2(3)
    @test QAtlas.fetch(w, ConformalWeights(); j=0) == 0 // 1
    @test QAtlas.fetch(w, ConformalWeights(); j=1 // 2) == 3 // 20    # (1/2)(3/2) / 5
    @test QAtlas.fetch(w, ConformalWeights(); j=1) == 2 // 5          # 1·2 / 5
    @test QAtlas.fetch(w, ConformalWeights(); j=3 // 2) == 3 // 4     # (3/2)(5/2) / 5 = 15/20
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=2)
end

@testset "WZWSU2: k=4 primary weights" begin
    w = WZWSU2(4)
    @test QAtlas.fetch(w, ConformalWeights(); j=2) == 1 // 1          # 2·3 / 6
    @test QAtlas.fetch(w, ConformalWeights(); j=1) == 1 // 3          # 1·2 / 6
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=5 // 2)
end

# ═══════════════════ Validation: half-integer / non-negative ═════════════════

@testset "WZWSU2: invalid j → DomainError" begin
    w = WZWSU2(4)
    # Negative j.
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=-1 // 2)
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=-1)
    # Non-half-integer (e.g. 1/3 is not a half-integer).
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=1 // 3)
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=2 // 3)
    # Float (not Integer / Rational).
    @test_throws DomainError QAtlas.fetch(w, ConformalWeights(); j=0.5)
end

# ═══════════════════ Type stability ═══════════════════════════════════════════

@testset "WZWSU2: ConformalWeights returns Rational{Int}" begin
    @test QAtlas.fetch(WZWSU2(1), ConformalWeights(); j=0) isa Rational{Int}
    @test QAtlas.fetch(WZWSU2(1), ConformalWeights(); j=1 // 2) isa Rational{Int}
    @test QAtlas.fetch(WZWSU2(7), ConformalWeights(); j=3) isa Rational{Int}
end
