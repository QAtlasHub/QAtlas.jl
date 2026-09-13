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

    # It is not a TFIM, so no TFIM method dispatches. What answers is not the
    # type system, though: AbstractQAtlas's top-level fallback catches every
    # unimplemented (model, quantity, bc) triple and names it. Pinning the type
    # said MethodError and got ErrorException, which is how that was found, so
    # pin the message instead.
    @test !(m isa TFIM)
    @test_throws "no fetch method for model=" fetch(m, MassGap(), Infinite())
end

@testset "Disordered :: the correlation decides which criterion applies" begin
    # Same chain, same family, same ν₀. Only the spatial correlation differs, and
    # the verdict flips. This is what makes the correlation an axis rather than a
    # label.
    fam = (; J=PowerLawDisorder(1.0), h=PowerLawDisorder(1.0))
    uncorr = Disordered(TFIM(), fam, Uncorrelated())
    fib = Disordered(TFIM(), fam, AperiodicSequence(-1.0))
    @test disorder_relevance(uncorr; d=1, d_euclidean=2) === :relevant   # ν₀ = 1 < 2/d
    @test disorder_relevance(fib; d=1, d_euclidean=2) === :irrelevant    # 1 > 1/(1−(−1))

    # ω = 1/2 is a random sequence, where Luck IS Harris in one dimension. The
    # wiring has to actually read corr.ω for this to mean anything, so check the
    # verdict MOVES with ω across Luck's own boundary at ν₀ = 1/(1−ω), i.e. ω = 0.
    rand_seq = Disordered(TFIM(), fam, AperiodicSequence(0.5))
    @test disorder_relevance(rand_seq; d=1, d_euclidean=2) ===
        disorder_relevance(uncorr; d=1, d_euclidean=2)
    for (ω, want) in ((-1.0, :irrelevant), (0.0, :marginal), (0.5, :relevant))
        m = Disordered(TFIM(), fam, AperiodicSequence(ω))
        @test disorder_relevance(m; d=1, d_euclidean=2) === want
    end

    # Correlated disorder reads the DISORDERED ν, which the atlas cannot supply.
    @test_throws "Pass `ν_dis`" disorder_relevance(
        Disordered(TFIM(), fam, PowerLawCorrelated(1.0)); d=1
    )
    # and each exponent is read by one criterion only: handing a route the other
    # one is refused rather than quietly dropped, which is how it used to go.
    @test_throws "does not read it" disorder_relevance(
        Disordered(TFIM(), fam, PowerLawCorrelated(1.0)); d=1, ν₀=1, ν_dis=2
    )
    @test_throws "only the Weinrib-Halperin criterion reads it" disorder_relevance(
        uncorr; d=1, d_euclidean=2, ν_dis=2
    )
    @test_throws "only the Weinrib-Halperin criterion reads it" disorder_relevance(
        fib; d=1, d_euclidean=2, ν_dis=2
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
    # `CriticalExponents`. There is no default, because `d + 1` is right for a
    # quantum chain and wrong for a classical model, and being wrong there is
    # silent: it answers, with the other model's number.
    @test_throws "no default" disorder_relevance(RandomTFIM(); d=1)
    @test disorder_relevance(RandomTFIM(); d=1, d_euclidean=2) === :relevant
    @test_throws ErrorException disorder_relevance(RandomTFIM(); d=1, d_euclidean=1)

    # Supplying ν₀ skips the lookup, so the Euclidean dimension stops mattering
    # and is not asked for.
    @test disorder_relevance(RandomTFIM(); d=1, ν₀=1) === :relevant
    @test disorder_relevance(RandomTFIM(); d=1, ν₀=3) === :irrelevant
    # Harris is marginal at ν₀ = 2/d exactly, so :marginal is reachable and not
    # a verdict the code can only ever name.
    @test disorder_relevance(RandomTFIM(); d=1, ν₀=2) === :marginal
    @test disorder_relevance(RandomTFIM(); d=2, ν₀=1) === :marginal
end

@testset "Disordered :: a missing atlas entry is named, not a MethodError" begin
    # The generic consumer works where the atlas knows the clean class and stops
    # with the specific gap where it does not. Both gaps are real today and the
    # messages tell them apart.
    no_class = Disordered(S1Heisenberg1D(); J=PowerLawDisorder(1.0))
    msg = try
        disorder_relevance(no_class; d=1, d_euclidean=2)
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("`UniversalityClass`", msg)
    @test occursin("S1Heisenberg1D", msg)
    # The atlas's own diagnosis is quoted, not replaced by a guess at what it was.
    @test occursin("The atlas said:", msg)

    # XXZ1D has a class, and that class's exponent set has no ν.
    no_nu = Disordered(XXZ1D(; Δ=0.5); Δ=PowerLawDisorder(1.0))
    msg2 = try
        disorder_relevance(no_nu; d=1, d_euclidean=2)
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("carries no", msg2)
    @test occursin("XY", msg2)                      # says which class
    # ...and does not tell the caller to fill the entry in, which is wrong advice
    # for a class whose transition has no power-law exponents to report.
    @test occursin("Berezinskii", msg2)
    @test !occursin("fill that entry in", msg2)
    @test msg != msg2                               # the two gaps differ
end
