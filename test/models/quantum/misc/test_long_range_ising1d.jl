# test/models/quantum/misc/test_long_range_ising1d.jl
#
# Phase 1 tests for LongRangeIsing1D (#293):
# only the α = Inf NN-TFIM limit is implemented (delegated to TFIM).
#
# Standalone:
#   julia --project=test test/models/quantum/misc/test_long_range_ising1d.jl

using QAtlas, Test

@testset "LongRangeIsing1D — α=Inf delegate to TFIM (Phase 1)" begin
    # Default: α=Inf, J=1, h=1 → QCP, gap = 0
    @test QAtlas.fetch(LongRangeIsing1D(), MassGap(), Infinite()) ≈ 0.0 atol = 1e-12
    # Paramagnetic h > J
    @test QAtlas.fetch(LongRangeIsing1D(; J=1.0, h=2.0), MassGap(), Infinite()) ≈ 2.0
    # Ferromagnetic h < J
    @test QAtlas.fetch(LongRangeIsing1D(; J=2.0, h=0.5), MassGap(), Infinite()) ≈ 3.0
    # Delegation invariant: matches TFIM directly
    Δ_lr = QAtlas.fetch(LongRangeIsing1D(; J=1.5, h=0.7), MassGap(), Infinite())
    Δ_tfim = QAtlas.fetch(TFIM(; J=1.5, h=0.7), MassGap(), Infinite())
    @test Δ_lr ≈ Δ_tfim
end

@testset "LongRangeIsing1D — finite α throws DomainError (Phase 2 deferral)" begin
    @test_throws DomainError QAtlas.fetch(LongRangeIsing1D(; α=2.0), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(LongRangeIsing1D(; α=10.0), MassGap(), Infinite())
    @test_throws DomainError QAtlas.fetch(LongRangeIsing1D(; α=0.5), MassGap(), Infinite())
end

@testset "LongRangeIsing1D — rejects J, h, α ≤ 0" begin
    @test_throws DomainError LongRangeIsing1D(; J=0.0)
    @test_throws DomainError LongRangeIsing1D(; h=-1.0)
    @test_throws DomainError LongRangeIsing1D(; α=0.0)
end
