# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Virasoro minimal models M(p, p_prime)
#
# Verifies the BPZ central-charge formula and the Kac-table conformal
# weights as exact `Rational{Int}` values, and the Kac symmetry
# h_{r,s} = h_{p_prime - r, p - s}.  Cross-checks the Ising special
# case M(4, 3) against `Universality(:Ising)`'s `c = 1//2` entry.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

# ═══════════════════ Construction & validation ════════════════════════════════

@testset "MinimalModel: construction validation" begin
    # Coprime requirement.
    @test_throws DomainError MinimalModel(6, 4)   # gcd 2
    @test_throws DomainError MinimalModel(9, 6)   # gcd 3
    @test_throws DomainError MinimalModel(10, 4)  # gcd 2
    # Order requirement: p > p_prime.
    @test_throws DomainError MinimalModel(3, 4)
    @test_throws DomainError MinimalModel(4, 4)   # also equal
    # Lower bound: p_prime >= 2.
    @test_throws DomainError MinimalModel(3, 1)
    @test_throws DomainError MinimalModel(2, 0)
    # Sanity: a few valid constructions return a struct with the right fields.
    m = MinimalModel(4, 3)
    @test m.p == 4 && m.p_prime == 3
end

# ═══════════════════ Central charge — exact rationals ═════════════════════════

@testset "MinimalModel: CentralCharge (exact)" begin
    @test QAtlas.fetch(MinimalModel(4, 3), CentralCharge()) == 1 // 2     # Ising
    @test QAtlas.fetch(MinimalModel(5, 4), CentralCharge()) == 7 // 10    # Tricritical Ising
    @test QAtlas.fetch(MinimalModel(6, 5), CentralCharge()) == 4 // 5     # 3-state Potts (chiral)
    @test QAtlas.fetch(MinimalModel(5, 2), CentralCharge()) == -22 // 5   # Yang–Lee
    @test QAtlas.fetch(MinimalModel(7, 6), CentralCharge()) == 6 // 7     # next unitary
    # Type stability: all of these are Rational{Int}.
    for m in
        (MinimalModel(4, 3), MinimalModel(5, 4), MinimalModel(6, 5), MinimalModel(5, 2))
        @test QAtlas.fetch(m, CentralCharge()) isa Rational{Int}
    end
end

# ═══════════════════ Cross-check: Ising = M(4, 3) ═════════════════════════════

@testset "MinimalModel: Ising cross-check (M(4,3) ↔ Universality(:Ising))" begin
    c_minimal = QAtlas.fetch(MinimalModel(4, 3), CentralCharge())
    c_ising = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).c
    @test c_minimal == c_ising
    @test c_minimal === 1 // 2
end

# ═══════════════════ Conformal weights — Ising primaries ══════════════════════

@testset "MinimalModel: Ising M(4,3) primary weights" begin
    m = MinimalModel(4, 3)
    # Ising primaries — identity, spin σ, energy ε.
    @test QAtlas.fetch(m, ConformalWeights(); r=1, s=1) == 0 // 1     # 1 (identity)
    @test QAtlas.fetch(m, ConformalWeights(); r=1, s=2) == 1 // 16    # σ
    @test QAtlas.fetch(m, ConformalWeights(); r=2, s=1) == 1 // 2     # ε
    # Type stability.
    @test QAtlas.fetch(m, ConformalWeights(); r=1, s=2) isa Rational{Int}
end

# ═══════════════════ Conformal weights — Tricritical Ising ════════════════════

@testset "MinimalModel: Tricritical Ising M(5,4) primary weights" begin
    m = MinimalModel(5, 4)
    # M(5, 4) primary weights computed from h_{r,s} = ((5r - 4s)^2 - 1) / 80
    # (Di Francesco–Mathieu–Sénéchal, Eq. 7.62 / Table 7.2).
    @test QAtlas.fetch(m, ConformalWeights(); r=1, s=1) == 0 // 1     # 1 (identity)
    @test QAtlas.fetch(m, ConformalWeights(); r=1, s=2) == 1 // 10    # ε
    @test QAtlas.fetch(m, ConformalWeights(); r=2, s=2) == 3 // 80    # σ
    @test QAtlas.fetch(m, ConformalWeights(); r=2, s=1) == 7 // 16    # σ′
    @test QAtlas.fetch(m, ConformalWeights(); r=1, s=3) == 3 // 5     # ε′
    @test QAtlas.fetch(m, ConformalWeights(); r=3, s=1) == 3 // 2     # ε′′
end

# ═══════════════════ Kac symmetry h_{r,s} = h_{p'-r, p-s} ═════════════════════

@testset "MinimalModel: Kac symmetry across many (p, p_prime)" begin
    for (p, q) in (
        (4, 3),    # Ising
        (5, 4),    # Tricritical Ising
        (6, 5),    # 3-state Potts (chiral)
        (7, 6),
        (8, 7),
        (5, 2),    # Yang–Lee (non-unitary)
        (7, 4),    # gcd = 1, p > p_prime
        (9, 4),
        (11, 6),
    )
        m = MinimalModel(p, q)
        for r in 1:(q - 1), s in 1:(p - 1)
            h_rs = QAtlas.fetch(m, ConformalWeights(); r=r, s=s)
            h_sym = QAtlas.fetch(m, ConformalWeights(); r=q - r, s=p - s)
            @test h_rs == h_sym
        end
    end
end

# ═══════════════════ Out-of-range (r, s) ══════════════════════════════════════

@testset "MinimalModel: out-of-range (r, s) → DomainError" begin
    m = MinimalModel(4, 3)
    # r ∈ 1..2, s ∈ 1..3 for M(4, 3).
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=0, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=3, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=1, s=0)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=1, s=4)
end

# ═══════════════════ PrimaryFields list ═══════════════════════════════════════

@testset "MinimalModel: PrimaryFields enumeration" begin
    # Ising M(4, 3): (p-1)(p_prime-1)/2 = 3*2/2 = 3 primaries (1, σ, ε).
    pf = QAtlas.fetch(MinimalModel(4, 3), PrimaryFields())
    @test length(pf) == 3
    hs = sort([x.h for x in pf])
    @test hs == [0 // 1, 1 // 16, 1 // 2]
    # Tricritical Ising M(5, 4): 4*3/2 = 6 primaries.
    pf2 = QAtlas.fetch(MinimalModel(5, 4), PrimaryFields())
    @test length(pf2) == 6
    @test all(x -> x.h isa Rational{Int}, pf2)
    # 3-state Potts (chiral) M(6, 5): 5*4/2 = 10 primaries.
    pf3 = QAtlas.fetch(MinimalModel(6, 5), PrimaryFields())
    @test length(pf3) == 10
end
