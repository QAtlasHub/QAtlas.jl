# The Griffiths dynamical exponent of the random transverse-field Ising chain.
#
# The implementation evaluates [(J/h)^{1/z}]_av = 1 (Iglói–Monthus Eq. (4.15),
# §4.1.3) in closed form on the model's disorder family. That evaluation is
# elementary and therefore exactly where an error would sit, so the tests below
# do not restate it: they put the returned z back into the DEFINING expectation
# by quadrature, and separately check the near-critical limit against the review's
# own closed form 1/z = 2|δ| (quoted with Eq. (4.51), §4.4.2).

using QAtlas, Test
using QAtlas: fetch, rtfim_delta
using QuadGK: quadgk

# E[(J/h)^{1/z}] for J = J·λ, h = h·μ with λ, μ ~ D⁻¹λ^{−1+1/D} on [0,1].
# Substituting λ = u^D makes u uniform, so each factor is a plain integral over
# [0,1] — computed here, not solved, so it cannot inherit the implementation's
# algebra.
#
# The inverse moment ∫₀¹ v^{−a} dv is only just integrable as a = D/z → 1, and
# deep in the disordered phase it gets there: z → D⁺, and at a = 0.978 the
# integrand overflows before quadgk can refine. A second substitution v = tᵐ
# turns it into ∫₀¹ t^{m(1−a)−1} dt, smooth once m(1−a) ≥ 1. `m` is chosen from
# `a`, which is an input to the check, so this changes the variable and not the
# question — verified against 1/(1−a) out to a = 0.999.
function _inverse_moment(a)
    m = max(2, ceil(Int, 2 / (1 - a)))
    val, _ = quadgk(t -> t^(m * (1 - a) - 1), 0, 1; rtol=1e-12)
    return m * val
end

function _griffiths_expectation(m::RandomTFIM, z)
    r = min(m.J, m.h) / max(m.J, m.h)
    pos, _ = quadgk(u -> u^(m.D / z), 0, 1; rtol=1e-12)
    return r^(1 / z) * pos * _inverse_moment(m.D / z)
end

@testset "RandomTFIM :: the returned z satisfies the defining expectation" begin
    for D in (0.5, 1.0, 2.0), hh in (1.25, 2.0, 5.0)
        m = RandomTFIM(; J=1.0, h=hh, D=D)
        z = fetch(m, DynamicalExponent(), Infinite())
        @test z > D                                   # the inverse moment must exist
        @test _griffiths_expectation(m, z) ≈ 1.0 rtol = 1e-8
    end
end

@testset "RandomTFIM :: z is dual under J ↔ h and grows toward criticality" begin
    # Interchanging h and J is the ordered-phase branch and must give the same z
    # (Iglói–Monthus, below Eq. (4.15)).
    for D in (0.7, 1.0), (a, b) in ((1.0, 3.0), (1.0, 1.4))
        @test fetch(RandomTFIM(; J=a, h=b, D=D), DynamicalExponent(), Infinite()) ≈
            fetch(RandomTFIM(; J=b, h=a, D=D), DynamicalExponent(), Infinite())
    end
    # Continuously varying, and divergent on approach — the property that makes it
    # a Griffiths exponent rather than a critical one.
    zs = [
        fetch(RandomTFIM(; J=1.0, h=hh, D=1.0), DynamicalExponent(), Infinite()) for
        hh in (4.0, 2.0, 1.5, 1.2, 1.05, 1.01)
    ]
    @test issorted(zs)                                # monotone in |δ|
    @test allunique(zs)                               # continuous, not a plateau
    @test zs[end] > 50                                # and running away
end

@testset "RandomTFIM :: near criticality 1/z → 2|δ|, with the next term" begin
    # The review states 1/z = 2|δ| near criticality. That is the LEADING term, so
    # rather than pick a threshold, expand: substituting x = 1/z into
    # exp(−2δD²x) = 1 − D²x² gives x = 2δ − 4δ³D², i.e. a RELATIVE error of
    # exactly 2δ²D². Asserting the coefficient pins the approach, not just that
    # it happens — a first-order-correct but second-order-wrong solver passes a
    # ratio band and fails this.
    for D in (0.5, 1.0, 2.0)
        for δ in (1e-2, 1e-3, 1e-4)
            m = RandomTFIM(; J=1.0, h=exp(2δ * D^2), D=D)
            @test rtfim_delta(m) ≈ δ rtol = 1e-12
            z = fetch(m, DynamicalExponent(), Infinite())
            @test 1 / z < 2δ                                  # approached from below
            @test abs(1 / z - 2δ) / (2δ) ≈ 2 * δ^2 * D^2 rtol = 1e-2
        end
    end
end

@testset "RandomTFIM :: criticality refuses z and supplies ψ instead" begin
    crit = RandomTFIM(; J=1.0, h=1.0, D=1.7)
    @test rtfim_delta(crit) == 0.0                    # whatever D — the cut-offs decide

    msg = try
        fetch(crit, DynamicalExponent(), Infinite())
        ""
    catch err
        err isa ErrorException ? sprint(showerror, err) : rethrow()
    end
    @test occursin("no finite dynamical exponent", msg)
    @test occursin("ActivatedExponent", msg)

    @test fetch(crit, ActivatedExponent(), Infinite()) === 1 // 2
    @test fetch(crit, UniversalityClass(), Infinite()) === Universality(:IsingSDRG)

    # ...and the exponents of that class are the ones already in the atlas.
    @test fetch(Universality(:IsingSDRG), CriticalExponents(); d=2).ψ ===
        fetch(crit, ActivatedExponent(), Infinite())
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
        @test occursin("δ = ", msg)                   # says how far off it is
    end
    # A finite z is what this chain does have.
    @test isfinite(fetch(off, DynamicalExponent(), Infinite()))
end

@testset "RandomTFIM :: the constructor rejects non-positive parameters" begin
    @test_throws ArgumentError RandomTFIM(; J=0.0, h=1.0, D=1.0)
    @test_throws ArgumentError RandomTFIM(; J=1.0, h=-1.0, D=1.0)
    @test_throws ArgumentError RandomTFIM(; J=1.0, h=1.0, D=0.0)
end
