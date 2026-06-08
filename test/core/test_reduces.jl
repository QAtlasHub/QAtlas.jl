# model -> model reductions (reduces / reduced_from) + C4 gap closure.

using QAtlas, Test
using QAtlas:
    TFIM,
    XXZ1D,
    Heisenberg1D,
    HeisenbergXYZ,
    MajumdarGhosh,
    J1J2Heisenberg1D,
    MixedFieldIsing1D,
    TricriticalPotts3,
    YangLee,
    KitaevHeisenberg,
    KitaevHoneycomb,
    S1XXZ1D,
    S1AnisotropicD1D,
    S1Heisenberg1D,
    MinimalModel,
    reductions,
    reduced_from

@testset "model -> model reductions (reduces)" begin
    @testset "reductions(model) lists the models a model reduces to" begin
        xyz = reductions(HeisenbergXYZ)
        @test any(r -> r.target === XXZ1D, xyz)
        @test all(r -> r.regime isa String && !isempty(r.regime), xyz)
        # a model can reduce to several targets (J1-J2 -> Heisenberg AND MajumdarGhosh)
        j1j2 = [r.target for r in reductions(J1J2Heisenberg1D)]
        @test Heisenberg1D in j1j2
        @test MajumdarGhosh in j1j2
    end

    @testset "reduced_from(model) lists the models reducing to it" begin
        into_tfim = [r.source for r in reduced_from(TFIM)]
        @test MixedFieldIsing1D in into_tfim          # h_z = 0
        @test !isempty(reduced_from(MajumdarGhosh))   # J1-J2 reduces to it at J2 = J1/2
        @test isempty(reduced_from(HeisenbergXYZ))    # nothing reduces TO the XYZ chain
    end

    @testset "reduces! appends a row" begin
        n = length(QAtlas.REDUCES)
        QAtlas.reduces!(HeisenbergXYZ, TFIM; regime="unit-test regime")
        @test length(QAtlas.REDUCES) == n + 1
        @test TFIM in [r.target for r in reductions(HeisenbergXYZ)]
    end
end

# The load-bearing assertion: every model→model delegation is now typed by a
# @reduces edge, so the C4 coherence check reports zero delegation gaps.
@testset "C4 delegation gaps closed by @reduces" begin
    fs = QAtlas.coherence_report()
    @test isempty(QAtlas.coherence_errors(fs))
    deleg_gaps = filter(g -> g.check === :delegation_target, QAtlas.coherence_gaps(fs))
    @test isempty(deleg_gaps)
end

# #658 — the C4/C5 gate and the `delegations` query must recognise descriptive
# `:<target>_delegation` method tags, not just the literal `:delegation`. Before
# this fix those rows silently bypassed C4, making the zero-gaps assertion above
# falsely green.
@testset "delegation gate recognizes variant method tags (#658)" begin
    @test QAtlas._is_delegation(:delegation)
    @test QAtlas._is_delegation(:kitaev_delegation)
    @test QAtlas._is_delegation(:minimal_model_delegation)
    @test QAtlas._is_delegation(:s1_heisenberg_delegation)
    @test !QAtlas._is_delegation(:analytic)
    @test !QAtlas._is_delegation(:closed_form)
end

@testset "variant-tagged model→model delegations are backed by @reduces (#658)" begin
    # each delegates via a :*_delegation tag and was previously invisible to C4;
    # the @reduces edge now types the delegation target.
    for (src, tgt) in (
        (TricriticalPotts3, MinimalModel),
        (YangLee, MinimalModel),
        (KitaevHeisenberg, KitaevHoneycomb),
        (S1XXZ1D, S1Heisenberg1D),
        (S1AnisotropicD1D, S1Heisenberg1D),
    )
        @test tgt in [r.target for r in reductions(src)]
    end
end
