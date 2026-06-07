# model <-> universality-class correspondence (realizes / realized_by).

using QAtlas, Test
using QAtlas: TFIM, XXZ1D, Heisenberg1D, realizations, realized_by

@testset "model <-> universality-class correspondence (realizes)" begin
    @testset "realizations(model) lists the classes a model realizes" begin
        tfim = realizations(TFIM)
        @test any(r -> r.class === :Ising, tfim)
        @test all(r -> r.regime isa String && !isempty(r.regime), tfim)
    end

    @testset "realized_by(class) lists member models" begin
        ising = [r.model for r in realized_by(:Ising)]
        @test TFIM in ising
        @test IsingSquare in ising            # quantum and 2D-classical both flow to Ising
        @test XXZ1D in [r.model for r in realized_by(:XY)]
        @test Heisenberg1D in [r.model for r in realized_by(:Heisenberg)]
        @test CurieWeissIsing in [r.model for r in realized_by(:MeanField)]
        @test TASEP in [r.model for r in realized_by(:KPZ)]
        @test isempty(realized_by(:NoSuchClass))
    end

    @testset "realizes! appends a row" begin
        n = length(QAtlas.REALIZES)
        QAtlas.realizes!(TFIM, :UnitTestClass; regime="unit-test regime")
        @test length(QAtlas.REALIZES) == n + 1
        @test only(realized_by(:UnitTestClass)).model === TFIM
    end
end
