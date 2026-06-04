# TFIM high-temperature free energy — the :approx worked example.
#
# Demonstrates the multi-definition framework for the :approx kind: the same
# (TFIM, FreeEnergy, Infinite) hub carries the exact closed form (canonical)
# AND a high-T expansion (scheme=:high_T), distinguishable by precision/domain.
# (verify_approx is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: TFIM, FreeEnergy, Infinite

@testset "TFIM FreeEnergy :approx — high-T expansion (scheme=:high_T)" begin
    m = TFIM(; J=1.0, h=0.5)
    f_ht(β) = -log(2) / β - (β / 2) * (m.J^2 + m.h^2)

    @testset "canonical (default) reproduces the exact closed form" begin
        f_exact = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.05)
        # bare fetch == scheme=:canonical, both = the existing exact row (goal 7)
        @test QAtlas.fetch(m, FreeEnergy(), Infinite(); scheme=:canonical, beta=0.05) ==
            f_exact
    end

    @testset "scheme=:high_T returns the expansion" begin
        @test QAtlas.fetch(m, FreeEnergy(), Infinite(); scheme=:high_T, beta=0.05) ≈
            f_ht(0.05) atol = 1e-12
        @test_throws ErrorException QAtlas.fetch(
            m, FreeEnergy(), Infinite(); scheme=:nonsense, beta=0.05
        )
    end

    @testset "agrees with exact in-domain (small β), breaks down out-of-domain (large β)" begin
        f_exact_lo = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.05)
        @test QAtlas.fetch(m, FreeEnergy(), Infinite(); scheme=:high_T, beta=0.05) ≈
            f_exact_lo rtol = 1e-3
        f_exact_hi = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=5.0)
        @test !isapprox(f_ht(5.0), f_exact_hi; rtol=1e-2)
    end

    @testset "definitions: exact (canonical) + high_T (:approx)" begin
        defs = QAtlas.definitions(m, FreeEnergy(), Infinite())
        @test length(defs) == 2
        canon = only(d for d in defs if d.canonical)
        @test canon.status === :exact
        ht = only(d for d in defs if d.scheme === :high_T)
        @test ht.status === :approx
        @test ht.valid_domain !== nothing && !isempty(ht.valid_domain)
        @test ht.references == ["Pfeuty1970"]
        v = QAtlas.validity(m, FreeEnergy(); scheme=:high_T)
        @test v.status === :approx
        @test v.error_order !== nothing
    end

    @testset "verify_approx card — in-domain agreement with the exact reference" begin
        β = 0.05
        s = verify_approx(
            m,
            FreeEnergy(),
            Infinite();
            route=:high_temperature,
            reference=QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β),
            agree_within=1e-3,
            valid_domain="betaJ << 1",
            error_order="O((betaJ)^3)",
            refs=["Pfeuty1970"],
            fetch_kw=(; scheme=:high_T, beta=β),
        )
        @test s ≈ f_ht(β) atol = 1e-12
    end
end
