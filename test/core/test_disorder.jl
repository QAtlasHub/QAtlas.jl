# Disorder decorates a model, and the decoration earns its place by having one
# generic consumer: `disorder_relevance` answers for ANY model, with no
# per-model code, because the correlation selects the criterion and the atlas
# supplies the clean exponent.

using QAtlas, Test
using QAtlas: fetch, clean_model, disorder, correlation

@testset "Disordered :: the decoration checks what it can" begin
    m = RandomTFIM(; J=1.0, h=2.0, D=1.0)
    @test m isa Disordered{TFIM}
    @test clean_model(m) == TFIM(; J=1.0, h=2.0)
    @test disorder(m, :J) == PowerLawDisorder(1.0)
    @test correlation(m) === Uncorrelated()

    # Any model with named couplings, which is the point of the decoration.
    @test Disordered(XXZ1D(; Δ=0.5); Δ=BinaryDisorder(0.4)) isa Disordered{XXZ1D}

    @test_throws ArgumentError Disordered(TFIM(); Jay=PowerLawDisorder(1.0))  # spelling
    @test_throws ArgumentError Disordered(TFIM(); J=1.0)                      # not a family
    @test_throws ArgumentError Disordered(TFIM())                             # nothing random
    @test_throws ArgumentError Disordered(TFIM(; J=0.0, h=1.0); J=PowerLawDisorder(1.0))
    # Not random is a different statement from random with no spread.
    @test_throws ArgumentError disorder(Disordered(TFIM(); J=PowerLawDisorder(1.0)), :h)

    # It is not a TFIM, so no clean method can dispatch on it. That is the type
    # system refusing, not a guard, which is why none is registered.
    @test !(m isa TFIM)
end

@testset "Disordered :: the correlation decides which criterion applies" begin
    # Same chain, same family, same ν₀. Only the spatial correlation differs, and
    # the verdict flips. This is what makes the correlation an axis rather than a
    # label.
    fam = (; J=PowerLawDisorder(1.0), h=PowerLawDisorder(1.0))
    uncorr = Disordered(TFIM(), fam, Uncorrelated())
    fib = Disordered(TFIM(), fam, Aperiodic(-1.0))
    @test disorder_relevance(uncorr; d=1) === :relevant      # Harris: ν₀ = 1 < 2/d
    @test disorder_relevance(fib; d=1) === :irrelevant       # Luck: 1 > 1/(1−(−1))

    # ω = 1/2 is a random sequence, where Luck IS Harris in one dimension.
    rand_seq = Disordered(TFIM(), fam, Aperiodic(0.5))
    @test disorder_relevance(rand_seq; d=1) === disorder_relevance(uncorr; d=1)

    # Correlated disorder reads the DISORDERED ν, which the atlas cannot supply.
    @test_throws ErrorException disorder_relevance(
        Disordered(TFIM(), fam, PowerLawCorrelated(1.0)); d=1
    )
    # Slower decay is more correlated, so more relevant: monotone in ρ through
    # the boundary at ρ = 2/ν_dis = 1.
    for (ρ, want) in ((0.5, :relevant), (1.0, :marginal), (2.0, :irrelevant))
        m = Disordered(TFIM(), fam, PowerLawCorrelated(ρ))
        @test disorder_relevance(m; d=1, ν_dis=2) === want
    end
end

@testset "Disordered :: the two dimensions are different numbers" begin
    # `d` is spatial and goes to the criterion; `d_euclidean` keys
    # `CriticalExponents` and defaults to d + 1. Passing one number to both is
    # the mistake the pair exists to prevent, and it fails loudly rather than
    # answering.
    @test disorder_relevance(RandomTFIM(); d=1) === :relevant
    @test disorder_relevance(RandomTFIM(); d=1, d_euclidean=2) === :relevant
    @test_throws ErrorException disorder_relevance(RandomTFIM(); d=1, d_euclidean=1)

    # Supplying ν₀ skips the lookup, so the Euclidean dimension stops mattering.
    @test disorder_relevance(RandomTFIM(); d=1, ν₀=1) === :relevant
    @test disorder_relevance(RandomTFIM(); d=1, ν₀=3) === :irrelevant
end

@testset "Disordered :: a missing atlas entry is named, not a MethodError" begin
    # The generic consumer works where the atlas knows the clean class and stops
    # with the specific gap where it does not. Both gaps are real today and the
    # messages tell them apart.
    no_class = Disordered(S1Heisenberg1D(); J=PowerLawDisorder(1.0))
    msg = try
        disorder_relevance(no_class; d=1)
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("no registered `UniversalityClass`", msg)

    # XXZ1D has a class, and that class's exponent set has no ν.
    no_nu = Disordered(XXZ1D(; Δ=0.5); Δ=PowerLawDisorder(1.0))
    msg2 = try
        disorder_relevance(no_nu; d=1)
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("has no `ν`", msg2)
    @test occursin("XY", msg2)                      # says which class
    @test msg != msg2                               # ...and the two gaps differ
end
