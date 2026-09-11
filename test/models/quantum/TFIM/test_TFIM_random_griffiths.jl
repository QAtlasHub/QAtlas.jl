# The Griffiths dynamical exponent of the random transverse-field Ising chain.
#
# One check carries the file: the returned z is put back into the DEFINING
# expectation, [(J/h)^{1/z}]_av = 1 (Iglói–Monthus Eq. (4.15), §4.1.3), by
# quadrature — so it cannot inherit the implementation's algebra. The rest test
# properties that expectation does not fix: the near-critical coefficient, the
# refusals, and the duality.

using QAtlas, Test
using QAtlas: fetch, rtfim_delta
using QuadGK: quadgk

# λ = u^D makes u uniform, so each factor of E[(J/h)^{1/z}] is an integral over
# [0,1].  ∫v^{−a}dv is only just integrable as a = D/z → 1, and deep in the
# disordered phase it reaches 0.978 and overflows; v = tᵐ with m from a makes it
# smooth (checked against 1/(1−a) out to a = 0.999).  That changes the variable,
# not the question.
function _inverse_moment(a)
    m = max(2, ceil(Int, 2 / (1 - a)))
    val, _ = quadgk(t -> t^(m * (1 - a) - 1), 0, 1; rtol=1e-12)
    return m * val
end

function _griffiths_expectation(m::RandomTFIM, z)
    pos, _ = quadgk(u -> u^(m.D / z), 0, 1; rtol=1e-12)
    return (min(m.J, m.h) / max(m.J, m.h))^(1 / z) * pos * _inverse_moment(m.D / z)
end

@testset "RandomTFIM :: the returned z satisfies the defining expectation" begin
    for D in (0.5, 1.0, 2.0), hh in (1.25, 2.0, 5.0)
        m = RandomTFIM(; J=1.0, h=hh, D=D)
        z = fetch(m, DynamicalExponent(), Infinite())
        # Not implied by the line below: for z < D the checker's own substitution
        # silently returns a finite wrong number instead of diverging.
        @test z > D
        @test _griffiths_expectation(m, z) ≈ 1.0 rtol = 1e-8
    end
end

@testset "RandomTFIM :: z varies continuously and diverges at criticality" begin
    zs = [
        fetch(RandomTFIM(; J=1.0, h=hh, D=1.0), DynamicalExponent(), Infinite()) for
        hh in (4.0, 2.0, 1.5, 1.2, 1.05, 1.01)
    ]
    @test issorted(zs) && allunique(zs)     # continuously varying, not a plateau
    @test zs[end] > 50                      # and running away — a Griffiths exponent
    # Duality: interchanging h and J is the ordered branch and gives the same z
    # (Iglói–Monthus, below Eq. (4.15)).  Guards writing J/h for min/max.
    @test fetch(RandomTFIM(; J=1.0, h=3.0, D=0.7), DynamicalExponent(), Infinite()) ≈
        fetch(RandomTFIM(; J=3.0, h=1.0, D=0.7), DynamicalExponent(), Infinite())
end

@testset "RandomTFIM :: near criticality 1/z → 2|δ| with the next term" begin
    # 1/z = 2|δ| is the LEADING term. Substituting x = 1/z into exp(−2δD²x) =
    # 1 − D²x² gives x = 2δ − 4δ³D², a relative error of exactly 2δ²D². Pinning
    # that coefficient catches a solver correct at first order and wrong at second,
    # which any ratio band would pass.
    for D in (0.5, 1.0, 2.0), δ in (1e-2, 1e-4)
        z = fetch(
            RandomTFIM(; J=1.0, h=exp(2δ * D^2), D=D), DynamicalExponent(), Infinite()
        )
        @test 1 / z < 2δ                                       # approached from below
        @test abs(1 / z - 2δ) / (2δ) ≈ 2 * δ^2 * D^2 rtol = 1e-2
    end
end

@testset "RandomTFIM :: criticality refuses z and supplies ψ instead" begin
    crit = RandomTFIM(; J=1.0, h=1.0, D=1.7)
    @test rtfim_delta(crit) == 0.0          # whatever D — the cut-offs alone decide
    @test rtfim_delta(RandomTFIM(; J=1.0, h=exp(2 * 1e-3), D=1.0)) ≈ 1e-3 rtol = 1e-12

    msg = try
        fetch(crit, DynamicalExponent(), Infinite())
        ""
    catch err
        err isa ErrorException ? sprint(showerror, err) : rethrow()
    end
    @test occursin("no finite dynamical exponent", msg)
    @test occursin("ActivatedExponent", msg)          # says what to ask for instead

    @test fetch(crit, UniversalityClass(), Infinite()) === Universality(:IsingSDRG)
    # ψ agrees with the exponent table of the class it flows to.
    @test fetch(crit, ActivatedExponent(), Infinite()) ===
        fetch(Universality(:IsingSDRG), CriticalExponents(); d=2).ψ ===
        1 // 2
end

@testset "RandomTFIM :: ψ and the class are refused off criticality" begin
    off = RandomTFIM(; J=1.0, h=2.0, D=1.0)
    for q in (ActivatedExponent(), UniversalityClass())
        msg = try
            fetch(off, q, Infinite())
            ""
        catch err
            err isa ErrorException ? sprint(showerror, err) : rethrow()
        end
        @test occursin("δ = ", msg)         # and says how far off it is
    end
    @test isfinite(fetch(off, DynamicalExponent(), Infinite()))

    @test_throws ArgumentError RandomTFIM(; J=0.0, h=1.0, D=1.0)
    @test_throws ArgumentError RandomTFIM(; J=1.0, h=1.0, D=0.0)
end
