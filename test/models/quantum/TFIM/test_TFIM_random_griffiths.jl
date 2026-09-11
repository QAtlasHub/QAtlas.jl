# ─────────────────────────────────────────────────────────────────────────────
# Griffiths dynamical exponent of the random transverse-field Ising chain.
#
# verify()-first.  The value pins are cards: z against a root found on the
# DEFINING expectation by quadrature (never touching the implementation's
# closed form), and against the review's published near-critical 1/z = 2|δ|.
# Raw @test is kept only where verify() cannot model the outcome — the two
# refusals, the structural inequalities, and the next-order coefficient, which
# is a statement about the ERROR of the card above rather than a value.
# ─────────────────────────────────────────────────────────────────────────────

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

# The independent route: bisect the defining condition itself.  E − 1 is +∞ as
# z → D⁺ (the inverse moment diverges) and negative as z → ∞ (where (ln r)/z
# beats D²/z²), so the bracket needs no input from the implementation.
function _z_by_quadrature(m::RandomTFIM)
    lo, hi = m.D * (1 + 1e-9), m.D * 2
    while _griffiths_expectation(m, hi) > 1
        hi *= 2
    end
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        (hi - lo) <= 1e-13 * mid && break
        _griffiths_expectation(m, mid) > 1 ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end

@testset "RandomTFIM :: Griffiths z" begin
    for D in (0.5, 1.0, 2.0), hh in (1.25, 2.0, 5.0)
        m = RandomTFIM(; J=1.0, h=hh, D=D)
        z = verify(
            m,
            DynamicalExponent(),
            Infinite();
            route=:sum_rule,
            independent=_z_by_quadrature(m),
            agree_within=1e-8,
            refs=[
                "Igloi-Monthus 2005 Eq. (4.15), §4.1.3: z is the positive root of [(J/h)^{1/z}]_av = 1",
            ],
            at=("D=$D", "h=$hh"),
        )
        # Not implied by the card: for z < D the checker's own substitution
        # returns a finite wrong number instead of diverging, so both sides
        # would agree on nonsense.
        @test z > D
    end
end

@testset "RandomTFIM :: near criticality, against the published 1/z = 2|δ|" begin
    # 2|δ| is the LEADING term, so the card's tolerance is the next one rather
    # than a guess: x = 2δ − 4δ³D² gives |z − 1/(2δ)| ≈ δD².
    for D in (0.5, 1.0, 2.0), δ in (1e-2, 1e-4)
        m = RandomTFIM(; J=1.0, h=exp(2δ * D^2), D=D)
        z = verify(
            m,
            DynamicalExponent(),
            Infinite();
            route=:literature_value,
            independent=1 / (2δ),
            agree_within=1.5 * δ * D^2,
            refs=[
                "Igloi-Monthus 2005, with Eq. (4.51), §4.4.2: 1/z = 2|δ| near criticality"
            ],
            at=("D=$D", "delta=$δ"),
        )
        # ...and the deviation the card tolerates is not slack: it is exactly
        # 2δ²D², which no tolerance can express. A solver right at first order
        # and wrong at second passes the card and fails here.
        @test 1 / z < 2δ
        @test abs(1 / z - 2δ) / (2δ) ≈ 2 * δ^2 * D^2 rtol = 1e-2
    end
end

@testset "RandomTFIM :: z varies continuously and diverges at criticality" begin
    # Structural, so raw: verify() pins values, not orderings.
    zs = [
        fetch(RandomTFIM(; J=1.0, h=hh, D=1.0), DynamicalExponent(), Infinite()) for
        hh in (4.0, 2.0, 1.5, 1.2, 1.05, 1.01)
    ]
    @test issorted(zs) && allunique(zs)     # continuously varying, not a plateau
    @test zs[end] > 50                      # and running away — a Griffiths exponent
    # Duality: interchanging h and J is the ordered branch and gives the same z
    # (Igloi-Monthus, below Eq. (4.15)).  Guards writing J/h for min/max.
    @test fetch(RandomTFIM(; J=1.0, h=3.0, D=0.7), DynamicalExponent(), Infinite()) ≈
        fetch(RandomTFIM(; J=3.0, h=1.0, D=0.7), DynamicalExponent(), Infinite())
end

@testset "RandomTFIM :: criticality refuses z and supplies ψ instead" begin
    # Raw @test_throws / message pins: verify() does not model error outcomes.
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
