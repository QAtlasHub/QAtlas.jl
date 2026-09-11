# The exact exponent table of the 1D random transverse-field Ising chain's
# infinite-randomness fixed point, and the one quantity it must refuse.
#
# Iglói–Monthus Table 1d. The values are pinned here as literature; what the
# tests below add is that the table is INTERNALLY consistent — the scaling
# relations β = ν·x_m, β_s = ν·x_m_s, ν_typ = ν(1−ψ) and φ = (d−x_m)/ψ each
# recompute one entry from others, so a typo in any single number fails.

using QAtlas, Test
using QAtlas: fetch

@testset "IsingSDRG :: exact IRFP exponents" begin
    e = fetch(Universality(:IsingSDRG), CriticalExponents(); d=2)

    @testset "the rational entries are exactly rational" begin
        # ν = 2, ψ = 1/2, ν_typ = 1, β_s = 1, x_m_s = 1/2 are exact, and stay
        # exact — a Float64 here would quietly lose the property later checks rely on.
        @test e.ν === 2 // 1
        @test e.ψ === 1 // 2
        @test e.ν_typ === 1 // 1
        @test e.β_s === 1 // 1
        @test e.x_m_s === 1 // 2
    end

    @testset "the table recomputes itself" begin
        # β = ν x_m and β_s = ν x_m_s — the review states both scaling relations
        # hold for this table.
        @test e.β ≈ e.ν * e.x_m rtol = 1e-15
        @test e.β_s == e.ν * e.x_m_s

        # ν_typ = ν(1−ψ): the typical correlation length, from the average one
        # and the activated exponent.
        @test e.ν_typ == e.ν * (1 - e.ψ)
        @test e.ν_typ < e.ν                      # strictly smaller, since ψ > 0

        # φ = (d − x_m)/ψ with the SPATIAL dimension d = 1 — the golden mean,
        # recomputed rather than quoted.
        @test e.φ ≈ (1 - e.x_m) / e.ψ rtol = 1e-14
        @test e.φ ≈ (1 + sqrt(5)) / 2 rtol = 1e-14
        # ...and with QAtlas's Euclidean d = 2 it is NOT the golden mean, which is
        # the trap the docstring warns about.
        @test !isapprox((2 - e.x_m) / e.ψ, (1 + sqrt(5)) / 2; rtol=1e-6)
    end

    @testset "the absent exponents stay absent" begin
        # α, γ, δ, η are not in Table 1d and cannot be filled in from the usual
        # scaling relations, which assume power-law dynamic scaling. Returning
        # them would be inventing them.
        for k in (:α, :γ, :δ, :η)
            @test !haskey(e, k)
        end
    end

    @testset "only d = 2 (1+1D), matching the CentralCharge method" begin
        @test_throws ErrorException fetch(
            Universality(:IsingSDRG), CriticalExponents(); d=3
        )
        # The Infinite-bc forwarding agrees with the 2-arg form.
        @test fetch(Universality(:IsingSDRG), CriticalExponents(), Infinite()) == e
    end
end

@testset "IsingSDRG :: ψ is fetchable and z is refused" begin
    @test fetch(Universality(:IsingSDRG), ActivatedExponent()) === 1 // 2
    @test fetch(Universality(:IsingSDRG), ActivatedExponent(), Infinite()) === 1 // 2

    # No finite z exists here. The refusal must DIAGNOSE, not merely fail:
    # a large number would read as a measurement.
    msg = try
        fetch(Universality(:IsingSDRG), DynamicalExponent())
        ""
    catch err
        err isa ErrorException ? sprint(showerror, err) : rethrow()
    end
    @test occursin("no finite dynamical exponent", msg)
    @test occursin("ActivatedExponent", msg)      # says what to ask for instead
    @test occursin("Griffiths", msg)              # ...and where a finite z does live

    # The clean class is unaffected: this refusal is about the fixed point, not
    # about the quantity being unsupported everywhere.
    @test fetch(Universality(:IsingSDRG), CentralCharge(); d=2) ≈ log(2) / 2
end
