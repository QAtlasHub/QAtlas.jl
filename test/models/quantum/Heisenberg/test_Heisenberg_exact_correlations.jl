using Test
using QAtlas

@testset "Heisenberg1D Exact Infinite Correlation Functions" begin
    model = Heisenberg1D()
    bc = Infinite()
    
    # Exact analytical values from Sato, Shiroishi, Takahashi (2005)
    exact_vals = [
        -0.14771572685331508,
        0.06067976995643609,
        -0.05024862725723622,
        0.03465277698273894,
        -0.030890366644598544
    ]
    
    for r in 1:5
        val_zz = fetch(model, ZZCorrelation{:static}(), bc; i=1, j=1+r, beta=Inf)
        @test isapprox(val_zz, exact_vals[r]; atol=1e-12)
        
        val_xx = fetch(model, XXCorrelation{:connected}(), bc; i=1, j=1+r, beta=Inf)
        @test isapprox(val_xx, exact_vals[r]; atol=1e-12)
    end
    
    # Check NotImplemented for r >= 6
    @test_throws ErrorException fetch(model, ZZCorrelation{:static}(), bc; i=1, j=7, beta=Inf)
end

@testset "S1Heisenberg1D Exact Correlation Functions Fallback" begin
    model = S1Heisenberg1D()
    @test_throws ErrorException fetch(model, ZZCorrelation{:static}(), Infinite(); i=1, j=2, beta=Inf)
end
