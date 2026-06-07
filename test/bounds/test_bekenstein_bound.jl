# Bekenstein entropy bound — a holographic bound in the bounds/ namespace
# (Bound{:Holographic}), extracted from the former Universality(:QuantumMechanics) dump.
# (verify_bound is in scope via the suite's global include of test/util/verify.jl.)

using QAtlas, Test
using QAtlas: Bound, BekensteinBound, Infinite

@testset "Bound(:Holographic) BekensteinBound — S ≤ 2πRE" begin
    holo = Bound(:Holographic)

    @test QAtlas.fetch(holo, BekensteinBound(), Infinite(); R=1.0, E=1.0) ≈ 2π atol = 1e-12
    @test QAtlas.fetch(holo, BekensteinBound(), Infinite(); R=2.0, E=3.0) ≈ 12π atol = 1e-12
    @test_throws ArgumentError QAtlas.fetch(
        holo, BekensteinBound(), Infinite(); R=-1.0, E=1.0
    )

    @testset "saturated by a black hole" begin
        R, E = 1.0, 1.0
        s = verify_bound(
            holo,
            BekensteinBound(),
            Infinite();
            route=:saturating_constant,
            measured=[2π * R * E],
            relation=:leq,
            saturating=true,
            refs=["Bekenstein1981"],
            fetch_kw=(; R=R, E=E),
        )
        @test s ≈ 2π atol = 1e-12
    end

    @testset "registered as an upper :bound" begin
        d = only(QAtlas.definitions(holo, BekensteinBound(), Infinite()))
        @test d.status === :bound
        @test d.direction === :upper
        @test d.canonical
    end
end
