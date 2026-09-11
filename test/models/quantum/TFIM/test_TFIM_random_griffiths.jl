# ─────────────────────────────────────────────────────────────────────────────
# Griffiths dynamical exponent of the random transverse-field Ising chain.
#
# verify()-first.  The value pins are cards: z against a root found on the
# DEFINING expectation by quadrature (never touching the implementation's
# closed form), and against the review's published near-critical 1/z = 2|δ|.
# Raw @test is kept only where verify() cannot model the outcome: the two
# refusals, the structural inequalities, and the next-order coefficient, which
# is a statement about the ERROR of the card above rather than a value.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using QAtlas: fetch, rtfim_delta
using QuadGK: quadgk

# A family that implements nothing, to check the contract error rather than a
# bare MethodError.
struct _PartialFamily <: DisorderFamily end

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
    pos, _ = quadgk(u -> u^(m.bonds.D / z), 0, 1; rtol=1e-12)
    return (m.J / m.h)^(1 / z) * pos * _inverse_moment(m.fields.D / z)
end

# The independent route: bisect the defining condition itself.  E − 1 is +∞ as
# z → D⁺ (the inverse moment diverges) and negative as z → ∞ (where (ln r)/z
# beats D²/z²), so the bracket needs no input from the implementation.
function _z_by_quadrature(m::RandomTFIM)
    # `m.fields.D` and NOT `moment_floor(m.fields)`: sharing that call with the
    # implementation makes the card blind to a bug in it. Measured, a factor-2
    # error there leaves impl and "independent" agreeing to 1e-9 while both sit
    # 65% off the true root.
    fl = m.fields.D
    lo, hi = fl * (1 + 1e-9), fl * 2
    steps = 0
    while _griffiths_expectation(m, hi) > 1
        hi *= 2
        # On the ordered branch the expectation approaches 1 from ABOVE and never
        # crosses, so this loop does not terminate without a cap.
        (steps += 1) > 200 && error("_z_by_quadrature: no bracket below z = $hi")
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
        # Not implied by the card: below the field family's floor the checker's
        # own substitution returns a finite wrong number instead of diverging,
        # so both sides would agree on nonsense.
        @test z > -1 / moment_floor(m.fields)
    end
end

@testset "RandomTFIM :: near criticality, against the published 1/z = 2|δ|" begin
    # 2|δ| is the LEADING term and is family-independent. The next one is NOT: it
    # comes from the fourth cumulant of ln λ, and for these two families it has
    # OPPOSITE signs. Measured (1/z−2δ)/(2δ) at δ = 1e-2 is −2.0e-4 for
    # PowerLaw(1) and +2.4e-5 for Binary(0.3). Asserting "approached from below"
    # for both, as this testset first did, would have been a power-law fact
    # written as a general one.
    for (f, coeff) in (
        (PowerLawDisorder(0.5), v -> -2 * v),        # x = 2δ − 4δ³D², v = D²
        (PowerLawDisorder(1.0), v -> -2 * v),
        (PowerLawDisorder(2.0), v -> -2 * v),
        (BinaryDisorder(0.3), v -> +2 * v / 3),
        (BinaryDisorder(0.6), v -> +2 * v / 3),
    )
        v = var_log(f)
        for δ in (1e-2, 1e-3)
            m = RandomTFIM(1.0, exp(2δ * v), f, f)   # δ = ln(h/J) / (2 var)
            @test rtfim_delta(m) ≈ δ rtol = 1e-12
            z = verify(
                m,
                DynamicalExponent(),
                Infinite();
                route=:literature_value,
                independent=1 / (2δ),
                agree_within=1.5 * δ * max(v, 1e-3),
                refs=[
                    "Igloi-Monthus 2005, with Eq. (4.51), §4.4.2: 1/z = 2|δ| near criticality",
                ],
                at=("family=$f", "delta=$δ"),
            )
            @test (1 / z - 2δ) / (2δ) ≈ coeff(v) * δ^2 rtol = 5e-2
        end
    end
end

@testset "RandomTFIM :: arbitrarily close to criticality, without a cut-off" begin
    # Solving in u = 1/z means no upper cap on z. At h = nextfloat(1.0) the true
    # root is 4.5e15 and the earlier z-space search gave up and returned Inf,
    # indistinguishable from the honest "exactly critical, no root" answer.
    m = RandomTFIM(1.0, nextfloat(1.0), PowerLawDisorder(1.0), PowerLawDisorder(1.0))
    δ = rtfim_delta(m)
    @test 0 < δ < 1e-15
    z = fetch(m, DynamicalExponent(), Infinite())
    @test isfinite(z)
    @test z ≈ 1 / (2δ) rtol = 1e-8          # ...and it is the paper's value
end

# BinaryDisorder makes the expectation a four-term sum, so this route shares no
# code with `log_moment` (no logs, no quadrature, no series) and no bracket with
# the implementation either.
function _binary_expectation(m::RandomTFIM, z)
    r = m.J / m.h
    return sum((r * λ / μ)^(1 / z) for λ in (1.0, m.bonds.κ), μ in (1.0, m.fields.κ)) / 4
end

function _z_by_enumeration(m::RandomTFIM)
    lo, hi = 1e-6, 1e8
    for _ in 1:300
        mid = 0.5 * (lo + hi)
        (hi - lo) <= 1e-13 * mid && break
        _binary_expectation(m, mid) > 1 ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end

@testset "RandomTFIM :: a second disorder family, by exact enumeration" begin
    # The condition is distribution-free; only its evaluation is per-family. A
    # family with a different STRUCTURE, bounded away from zero so every moment
    # exists, checks that the machinery is general and not the power law in
    # disguise.
    for (κJ, κh, hh) in ((0.3, 0.3, 2.0), (0.3, 0.3, 1.1), (0.1, 0.1, 5.0))
        m = RandomTFIM(1.0, hh, BinaryDisorder(κJ), BinaryDisorder(κh))
        verify(
            m,
            DynamicalExponent(),
            Infinite();
            route=:sum_rule,
            independent=_z_by_enumeration(m),
            agree_within=1e-8,
            refs=["Igloi-Monthus 2005 Eq. (4.15), §4.1.3, on two-valued couplings"],
            at=("kappa=$κJ", "h=$hh"),
        )
    end
    # No z floor here, unlike the power law: `z > D` is a property of that family
    # and not of the condition.
    @test moment_floor(BinaryDisorder(0.3)) == -Inf
    @test moment_floor(PowerLawDisorder(2.0)) == -0.5
end

@testset "RandomTFIM :: the family interface is consistent and enforced" begin
    # `log_moment(f, s) = ln E[λ^s]` is the cumulant generating function of ln λ, so
    # `mean_log` and `var_log` are its first two derivatives at s = 0. They are
    # supplied in closed form for precision, since `rtfim_delta` is checked to
    # 1e-12 and a finite difference is nowhere near that, but nothing bound the three
    # together, so an algebra slip in one would not have shown up anywhere.
    #
    # ABSOLUTE tolerance, not relative: `var_log` → 0 as κ → 1 (no disorder), so a
    # relative bound is unbounded there. Measured, the relative error on
    # BinaryDisorder(0.99) is 2e-4 while the absolute one is 5e-9. At h = 1e-4 the
    # worst absolute error over these families is 4.2e-8, so 1e-6 is achievable
    # with room, and still catches any real algebra slip, which is O(1).
    for f in (
        PowerLawDisorder(1.7),
        PowerLawDisorder(0.4),
        BinaryDisorder(0.35),
        BinaryDisorder(0.9),
        BinaryDisorder(0.99),            # var_log ≈ 2.5e-5, the hard one
    )
        h = 1e-4
        d1 = (log_moment(f, h) - log_moment(f, -h)) / (2h)
        d2 = (log_moment(f, h) - 2 * log_moment(f, 0.0) + log_moment(f, -h)) / h^2
        @test mean_log(f) ≈ d1 atol = 1e-6
        @test var_log(f) ≈ d2 atol = 1e-6
    end

    # A family implementing only some of the four must SAY so. Three of the four
    # fetch routes never reach `moment_floor`, so without this a partial family
    # passes them and fails later with a bare MethodError.
    msg = try
        moment_floor(_PartialFamily())
        ""
    catch err
        err isa ErrorException ? sprint(showerror, err) : rethrow()
    end
    @test occursin("DisorderFamily", msg)
    @test occursin("all four", msg)
end

@testset "RandomTFIM :: bonds and fields may be different families" begin
    # The two type parameters are only worth having if the shapes can differ, not
    # merely the parameters. Every other test here varies κ or D within one family.
    m = RandomTFIM(1.0, 2.0, PowerLawDisorder(1.0), BinaryDisorder(0.3))
    @test m isa RandomTFIM{PowerLawDisorder,BinaryDisorder}
    z = fetch(m, DynamicalExponent(), Infinite())
    @test isfinite(z) && z > 0
    # Mixed E[(J/h)^{1/z}]: power-law bonds by quadrature, two-valued fields by sum.
    pos, _ = quadgk(u -> u^(m.bonds.D / z), 0, 1; rtol=1e-12)
    neg = (1.0^(-1 / z) + m.fields.κ^(-1 / z)) / 2
    @test (m.J / m.h)^(1 / z) * pos * neg ≈ 1.0 rtol = 1e-8
end

@testset "RandomTFIM :: no rare regions means no exponent, said so" begin
    # For two-valued couplings a root exists only while J·max λ > h·min μ: the
    # strongest bond must beat the weakest field, or nothing can be locally ordered
    # and there is no Griffiths phase. The sweep above sits ON that boundary:
    # (κ, h) = (0.1, 5.0) has J/h = 0.2 > κ and works, and κ = 0.2 at the same h
    # does not. Before the bracket was checked, both returned 1.00003e-8, which is
    # the bisection's own starting point, dressed as an exponent.
    @test fetch(
        RandomTFIM(1.0, 5.0, BinaryDisorder(0.1), BinaryDisorder(0.1)),
        DynamicalExponent(),
        Infinite(),
    ) ≈ 0.507823054 rtol = 1e-6
    for κ in (0.2, 0.3, 0.9)
        msg = try
            fetch(
                RandomTFIM(1.0, 5.0, BinaryDisorder(κ), BinaryDisorder(κ)),
                DynamicalExponent(),
                Infinite(),
            )
            ""
        catch err
            err isa ErrorException ? sprint(showerror, err) : rethrow()
        end
        @test occursin("has no root", msg)
        @test occursin("rare regions", msg)   # ...and says why, not just that
    end

    # A deterministic coupling is allowed by `BinaryDisorder` (0 < κ ≤ 1) but makes
    # the spread of ln λ zero, so δ is 0/0 or finite/0. `iszero` sees neither.
    for m in (
        RandomTFIM(1.0, 1.0, BinaryDisorder(1.0), BinaryDisorder(1.0)),
        RandomTFIM(1.0, 2.0, BinaryDisorder(1.0), BinaryDisorder(1.0)),
    )
        @test !isfinite(rtfim_delta(m))
        msg = try
            fetch(m, DynamicalExponent(), Infinite())
            ""
        catch err
            err isa ErrorException ? sprint(showerror, err) : rethrow()
        end
        @test occursin("no disorder", msg)
    end
end

@testset "RandomTFIM :: duality swaps the whole problem, not the two scales" begin
    # With unequal bond and field families, min(J,h)/max(J,h) is NOT the dual:
    # it leaves the ordered side without a root. The swap is (J, bonds) ↔
    # (h, fields) (Igloi-Monthus, below Eq. (4.15)).
    ordered = RandomTFIM(1.0, 1.0, BinaryDisorder(0.5), BinaryDisorder(0.2))
    dual = RandomTFIM(1.0, 1.0, BinaryDisorder(0.2), BinaryDisorder(0.5))
    @test rtfim_delta(ordered) ≈ -rtfim_delta(dual)
    @test rtfim_delta(ordered) < 0                      # the ordered branch
    @test fetch(ordered, DynamicalExponent(), Infinite()) ≈
        fetch(dual, DynamicalExponent(), Infinite())

    # Self-consistency alone is weak: a bonds/fields swap inside the residual makes
    # BOTH sides wrong identically. So check the ORDERED branch against ground
    # truth, by quadrature on the dual problem the implementation should be solving.
    ord = RandomTFIM(3.0, 1.0, PowerLawDisorder(0.7), PowerLawDisorder(0.7))
    @test rtfim_delta(ord) < 0
    @test fetch(ord, DynamicalExponent(), Infinite()) ≈ _z_by_quadrature(
        RandomTFIM(1.0, 3.0, PowerLawDisorder(0.7), PowerLawDisorder(0.7))
    ) rtol = 1e-8
    # ...and criticality is [ln J] = [ln h], which here is NOT J == h.
    @test !iszero(
        rtfim_delta(RandomTFIM(1.0, 1.0, BinaryDisorder(0.5), BinaryDisorder(0.2)))
    )
    @test iszero(
        rtfim_delta(RandomTFIM(1.0, 1.0, BinaryDisorder(0.5), BinaryDisorder(0.5)))
    )
end

@testset "RandomTFIM :: z varies continuously and diverges at criticality" begin
    # Structural, so raw: verify() pins values, not orderings.
    zs = [
        fetch(RandomTFIM(; J=1.0, h=hh, D=1.0), DynamicalExponent(), Infinite()) for
        hh in (4.0, 2.0, 1.5, 1.2, 1.05, 1.01)
    ]
    @test issorted(zs) && allunique(zs)     # continuously varying, not a plateau
    @test zs[end] > 50                      # and running away, a Griffiths exponent
    # Duality: interchanging h and J is the ordered branch and gives the same z
    # (Igloi-Monthus, below Eq. (4.15)).  Guards writing J/h for min/max.
    @test fetch(RandomTFIM(; J=1.0, h=3.0, D=0.7), DynamicalExponent(), Infinite()) ≈
        fetch(RandomTFIM(; J=3.0, h=1.0, D=0.7), DynamicalExponent(), Infinite())
end

@testset "RandomTFIM :: criticality refuses z and supplies ψ instead" begin
    # Raw @test_throws / message pins: verify() does not model error outcomes.
    crit = RandomTFIM(; J=1.0, h=1.0, D=1.7)
    @test rtfim_delta(crit) == 0.0          # whatever D; the cut-offs alone decide
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
