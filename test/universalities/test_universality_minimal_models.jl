# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Virasoro minimal models M(p, p_prime)
#
# Migrated from pure-legacy @test to verify()-first (PR #449 phase B,
# zero-legacy end-state). BPZ central charges and Kac-table conformal
# weights become verify() cards (:second_closed_form). Construction errors,
# type-stability traits, Kac-symmetry sweep, and PrimaryFields enumeration
# stay raw @test (verify() can't represent error-throwing dispatches,
# type-isa checks, or multi-element invariants).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "MinimalModel: construction validation" begin
    @test_throws DomainError MinimalModel(6, 4)
    @test_throws DomainError MinimalModel(9, 6)
    @test_throws DomainError MinimalModel(10, 4)
    @test_throws DomainError MinimalModel(3, 4)
    @test_throws DomainError MinimalModel(4, 4)
    @test_throws DomainError MinimalModel(3, 1)
    @test_throws DomainError MinimalModel(2, 0)
    m = MinimalModel(4, 3)
    @test m.p == 4 && m.p_prime == 3
end

@testset "MinimalModel: CentralCharge (BPZ exact rationals)" begin
    for (p, q, c_lit, label) in (
        (4, 3, 1 // 2, "Ising"),
        (5, 4, 7 // 10, "Tricritical Ising"),
        (6, 5, 4 // 5, "3-state Potts (chiral)"),
        (5, 2, -22 // 5, "Yang-Lee (non-unitary)"),
        (7, 6, 6 // 7, "next unitary"),
    )
        verify(
            MinimalModel(p, q),
            CentralCharge(),
            Infinite();
            route=:second_closed_form,
            independent=c_lit,
            agree_within=0,
            at=["M($(p),$(q)) = $(label)"],
            refs=[
                "BPZ 1984 (Belavin-Polyakov-Zamolodchikov): c = 1 - 6(p-p_prime)^2/(p*p_prime); M($(p),$(q)) ⇒ c = $(c_lit)",
            ],
        )
    end
    for m in
        (MinimalModel(4, 3), MinimalModel(5, 4), MinimalModel(6, 5), MinimalModel(5, 2))
        @test QAtlas.fetch(m, CentralCharge()) isa Rational{Int}
    end
end

@testset "MinimalModel: Ising cross-check (M(4,3) ↔ Universality(:Ising))" begin
    c_minimal = QAtlas.fetch(MinimalModel(4, 3), CentralCharge())
    c_ising = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).c
    @test c_minimal == c_ising
    @test c_minimal === 1 // 2
end

@testset "MinimalModel: Ising M(4,3) primary weights" begin
    verify(
        MinimalModel(4, 3),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=0 // 1,
        agree_within=0,
        at=["(r=1,s=1) identity"],
        refs=["BPZ 1984 Kac table: M(4,3) (1,1) ⇒ 0"],
        fetch_kw=(; r=1, s=1),
    )
    verify(
        MinimalModel(4, 3),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=1 // 16,
        agree_within=0,
        at=["(r=1,s=2) σ spin"],
        refs=["BPZ Kac table: M(4,3) (1,2) ⇒ h_σ = 1/16"],
        fetch_kw=(; r=1, s=2),
    )
    verify(
        MinimalModel(4, 3),
        ConformalWeights(),
        Infinite();
        route=:second_closed_form,
        independent=1 // 2,
        agree_within=0,
        at=["(r=2,s=1) ε energy"],
        refs=["BPZ Kac table: M(4,3) (2,1) ⇒ h_ε = 1/2"],
        fetch_kw=(; r=2, s=1),
    )
    @test QAtlas.fetch(MinimalModel(4, 3), ConformalWeights(); r=1, s=2) isa Rational{Int}
end

@testset "MinimalModel: Tricritical Ising M(5,4) primary weights" begin
    for (r, s, h_lit, label) in (
        (1, 1, 0 // 1, "1 identity"),
        (1, 2, 1 // 10, "ε"),
        (2, 2, 3 // 80, "σ"),
        (2, 1, 7 // 16, "σ′"),
        (1, 3, 3 // 5, "ε′"),
        (3, 1, 3 // 2, "ε′′"),
    )
        verify(
            MinimalModel(5, 4),
            ConformalWeights(),
            Infinite();
            route=:second_closed_form,
            independent=h_lit,
            agree_within=0,
            at=["(r=$(r),s=$(s)) $(label)"],
            refs=[
                "DFMS 1997 Eq. 7.62 / Table 7.2: M(5,4) Tricritical Ising primary weights h_{r,s} = ((5r - 4s)^2 - 1)/80",
            ],
            fetch_kw=(; r=r, s=s),
        )
    end
end

@testset "MinimalModel: Kac symmetry across many (p, p_prime)" begin
    for (p, q) in ((4, 3), (5, 4), (6, 5), (7, 6), (8, 7), (5, 2), (7, 4), (9, 4), (11, 6))
        m = MinimalModel(p, q)
        for r in 1:(q - 1), s in 1:(p - 1)
            h_rs = QAtlas.fetch(m, ConformalWeights(); r=r, s=s)
            h_sym = QAtlas.fetch(m, ConformalWeights(); r=q - r, s=p - s)
            @test h_rs == h_sym
        end
    end
end

@testset "MinimalModel: out-of-range (r, s) → DomainError" begin
    m = MinimalModel(4, 3)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=0, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=3, s=1)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=1, s=0)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(); r=1, s=4)
end

@testset "MinimalModel: PrimaryFields enumeration" begin
    pf = QAtlas.fetch(MinimalModel(4, 3), PrimaryFields())
    @test length(pf) == 3
    hs = sort([x.h for x in pf])
    @test hs == [0 // 1, 1 // 16, 1 // 2]
    pf2 = QAtlas.fetch(MinimalModel(5, 4), PrimaryFields())
    @test length(pf2) == 6
    @test all(x -> x.h isa Rational{Int}, pf2)
    pf3 = QAtlas.fetch(MinimalModel(6, 5), PrimaryFields())
    @test length(pf3) == 10
end
