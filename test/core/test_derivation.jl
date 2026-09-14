# test/core/test_derivation.jl — the scaling-plane generator itself
# (src/core/derivation.jl), as opposed to the checks it emits
# (test/generated/test_derivation_checks.jl).
#
# The load-bearing testset here is the mutation control.  A cross-check whose
# tolerance is derived from the table's own quoted errors can be made to pass
# everything by a tolerance that is too wide, and nothing about a green run says
# which.  So: move one exponent and count the routes that go red.

using QAtlas, Test
using QAtlas:
    _scaling_data,
    _route_sigma,
    _sigma_outcome,
    _roundoff_floor,
    exponent_sweep!,
    refuse_exponents!,
    exponent_hubs,
    check_derivation_coverage,
    derivation_reach,
    EXPONENT_SWEEPS
using AbstractQAtlas: consistency_report

# The 2D Ising table, exact; and the 3D percolation table with its quoted errors.
const _DV_ISING2D = (α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4)
const _DV_PERC3D = (
    α=-0.625,
    α_err=0.003,
    β=0.4181,
    β_err=0.0008,
    γ=1.793,
    γ_err=0.003,
    δ=5.29,
    δ_err=0.06,
    ν=0.87619,
    ν_err=0.00012,
    η=-0.04590,
    η_err=0.00030,
)

# How many (target, route) pairs of `table` land outside k sigma, with sigma the
# quoted error of the target and the route's inputs propagated through it.  The
# generator's verdict, reproduced here on data the atlas need not hold.
function _dv_red(table::NamedTuple, dim; k=3.0)
    data, errs = _scaling_data(table, dim)
    n_red = 0
    n_all = 0
    for row in consistency_report(data; domain=:scaling, atol=0, rtol=0)
        for (step, value) in zip(row.steps, row.values)
            n_all += 1
            out = _sigma_outcome(
                row.held_out,
                value,
                getfield(errs, row.target),
                _route_sigma(step, data, errs, value);
                k=k,
            )
            out.status === :fail && (n_red += 1)
        end
    end
    return (n_red, n_all)
end

@testset "an exact table agrees on every route, to round-off" begin
    red, all = _dv_red(_DV_ISING2D, 2 // 1)
    @test all == 12                 # 7 held-out variables, 12 routes back
    @test red == 0
end

@testset "mutation control: the routes have teeth on this data" begin
    # Exact table: a 5% move in any exponent must go red somewhere.  Reported per
    # exponent rather than summed, because a sum can be carried by one sensitive
    # exponent while another is unconstrained, and that is exactly the blind spot
    # worth naming.
    for sym in (:α, :β, :γ, :δ, :ν, :η)
        base = getfield(_DV_ISING2D, sym)
        moved = iszero(base) ? 1 // 20 : base * 21 // 20
        bad = merge(_DV_ISING2D, NamedTuple{(sym,)}((moved,)))
        red, all = _dv_red(bad, 2 // 1)
        @test red ≥ 3
        @test red < all              # a single wrong exponent is not every route
    end
end

@testset "a quoted-error table is judged on its own errors, not a guessed rtol" begin
    red, all = _dv_red(_DV_PERC3D, 3.0)
    @test all == 12
    @test red == 0

    # The defect this plane was built on: η of the wrong sign and magnitude.  It
    # is the ONLY thing moved, and the Fisher routes deny it.
    wrong = merge(_DV_PERC3D, (η=0.46, η_err=0.08))
    red_w, _ = _dv_red(wrong, 3.0)
    @test red_w ≥ 3

    # ...and a value inside the quoted errors is not a failure, however far it is
    # in relative terms: a 1-sigma move of a small exponent must stay green.
    nudged = merge(_DV_PERC3D, (η=_DV_PERC3D.η + _DV_PERC3D.η_err,))
    red_n, _ = _dv_red(nudged, 3.0)
    @test red_n == 0
end

@testset "sigma propagates through the route, not just the target" begin
    data, errs = _scaling_data(_DV_PERC3D, 3.0)
    rows = consistency_report(data; domain=:scaling, atol=0, rtol=0)
    row = only(filter(r -> r.target === :η, rows))
    step, value = only(row.steps), only(row.values)
    @test nameof(typeof(step.relation)) === :Fisher

    # Fisher solved for η is 2 − γ/ν, so the route's sigma is γ_err/ν and
    # γ·ν_err/ν² in quadrature.  Hand-computed, because a propagation that
    # silently returned zero would make every route look infinitely precise.
    expected = hypot(
        _DV_PERC3D.γ_err / _DV_PERC3D.ν, _DV_PERC3D.γ * _DV_PERC3D.ν_err / _DV_PERC3D.ν^2
    )
    # `rtol` is loose because the propagation is a finite difference AT the quoted
    # error and Fisher is not affine in ν; what is asserted is the scale, not a
    # digit count.
    @test _route_sigma(step, data, errs, value) ≈ expected rtol = 1e-4
    @test _route_sigma(step, data, errs, value) > 10 * _DV_PERC3D.η_err
end

@testset "an inexact value with no quoted error is refused, not judged" begin
    # Drop one `_err` and the table can no longer state how precise it is.
    stripped = NamedTuple{filter(!=(:γ_err), keys(_DV_PERC3D))}(
        Tuple(_DV_PERC3D[k] for k in keys(_DV_PERC3D) if k !== :γ_err)
    )
    out = _scaling_data(stripped, 3.0)
    @test out isa String
    @test occursin("γ_err", out)

    # Too few exponents to close is likewise a refusal with a reason.
    @test _scaling_data((; η=1 // 4), 2 // 1) isa String
    @test occursin("too few to close", _scaling_data((; η=1 // 4), 2 // 1))
end

@testset "the round-off floor is machine precision, never a physics tolerance" begin
    @test _roundoff_floor(1.0, 1.0) < 1e-14
    @test _roundoff_floor(1e6, 1e6) < 1e-8
    # An exact pair passes with zero sigma; a pair one part in 10^6 apart does not.
    @test _sigma_outcome(1.75, 1.75, 0, 0; k=3).status === :pass
    @test _sigma_outcome(1.75, 1.750002, 0, 0; k=3).status === :fail
    # `rel_err` carries the deviation in sigma, which is what a failing row needs.
    @test _sigma_outcome(1.0, 1.3, 0.1, 0.0; k=3).rel_err ≈ 3.0
end

@testset "every exponent hub is declared, and the guard can see them" begin
    # Positive control first: a guard reading an empty hub list would report no
    # gaps for the wrong reason, and `isempty(findings)` cannot tell the two
    # apart.  Fifteen types carry a CriticalExponents method: the ten classes and
    # five aliases below, plus the four models that delegate to them.
    hubs = exponent_hubs()
    @test length(hubs) ≥ 15
    for M in (
        Universality{:Ising},
        Universality{:Percolation},
        Universality{:Potts3},
        Universality{:Potts4},
        Universality{:XY},
        Universality{:Heisenberg},
        Universality{:MeanField},
        Universality{:IsingSDRG},
        QAtlas.Ising2D,
        QAtlas.MeanField,
        QAtlas.KPZ1D,
        # ...and the delegating models, which the coverage guard must SEE in order
        # to exempt them; an unseen hub and an exempted one both report no gap.
        QAtlas.TFIM,
        QAtlas.IsingSquare,
        QAtlas.IsingTriangular,
        QAtlas.CurieWeissIsing,
    )
        @test M in hubs
    end
    @test isempty(check_derivation_coverage())
end

@testset "a declaration that would read as coverage is refused" begin
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Ising}, Infinite; sweep=(d=[2],)
    )                                                  # duplicate (model, bc, sweep)
    @test_throws ArgumentError exponent_sweep!(Universality{:Potts3}, PBC; k=0)
    @test_throws ArgumentError exponent_sweep!(Universality{:Potts3}, PBC; sweep=(;))
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Potts3}, PBC; dimension=:nope, sweep=(d=[2],)
    )
    # A non-numeric swept `d` reaches the relations without passing the `isa Real`
    # test the fetched exponents get; caught here, it names the declaration.
    @test_throws ArgumentError exponent_sweep!(Universality{:Potts3}, PBC; sweep=(d=["2"],))
    # dimension=nothing suppresses the hyperscaling routes; unexplained, that is
    # a silent exclusion wearing the shape of a pass.
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Potts3}, PBC; dimension=nothing
    )
    @test_throws ArgumentError refuse_exponents!(Universality{:Potts3}, PBC; reason="")
    # A refusal is keyed on (model, bc) alone: its check id carries no sweep point,
    # so a second one would pass a sweep-aware test and then collide on that id.
    @test_throws ArgumentError refuse_exponents!(KPZ1D, Infinite; reason="again")
    # ...and nothing above left a row behind.
    @test !any(s -> s.bc === PBC, EXPONENT_SWEEPS)
    @test count(s -> s.model === KPZ1D, EXPONENT_SWEEPS) == 1
end

# `derived_from` is a second copy of a provenance header, and it fails in the
# direction that looks fine: a name matching nothing suppresses nothing, the
# route runs, and it PASSES, because the shipped value is what that relation
# predicts. Both halves are validated at declaration.
@testset "derived_from cannot name an exponent or a relation that does not exist" begin
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Potts3}, PBC; sweep=(d=[2],), derived_from=[:delta => [:Widom]]
    )
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Potts3}, PBC; sweep=(d=[2],), derived_from=[:δ => [:widom]]
    )
    @test !any(s -> s.bc === PBC, EXPONENT_SWEEPS)

    # Positive control: the spelling the registry actually uses is accepted, so
    # the two refusals above are rejecting the typo and not the shape.
    known = QAtlas._scaling_relation_names()
    for r in (:Widom, :Fisher, :Rushbrooke, :Josephson)
        @test r in known
    end
end

# The emptiness test that reads naturally cannot fire: a point that cannot be
# prepared still emits a refusal check. What CAN happen is a declaration all of
# whose points were refused, and that is what is asked.
@testset "a declaration that judges nothing is reported, and the guard can fire" begin
    for s in EXPONENT_SWEEPS
        s isa QAtlas.SweptExponents || continue
        @test !QAtlas._emits_only_skips(s)
    end
    # Positive control, built directly so it never enters the store: the BKT point
    # returns η alone, which is too few to close, so every point is a skip.
    all_skip = QAtlas.SweptExponents(
        Universality{:XY},
        Infinite,
        (d=[2],),
        :sweep,
        Pair{Symbol,Vector{Symbol}}[],
        3.0,
        "",
        String[],
    )
    @test QAtlas._emits_only_skips(all_skip)

    # A hub whose fetch does not work is NOT that case: it is a config bug and is
    # reported `:error`, which fails the suite on its own. A coverage finding on
    # top would be a second report of one thing, and filing it as a skip would be
    # the silence this whole split exists to prevent.
    dead = QAtlas.SweptExponents(
        QAtlas.Heisenberg1D,
        Infinite,
        (;),
        2,
        Pair{Symbol,Vector{Symbol}}[],
        3.0,
        "",
        String[],
    )
    checks = QAtlas._scaling_checks(dead)
    @test !isempty(checks)
    @test all(c -> QAtlas.run_generated_check(c).status === :error, checks)
    @test !QAtlas._emits_only_skips(dead)
end

@testset "a non-finite sigma is refused, never used as a tolerance" begin
    # The tolerance is built FROM the sigma, so an infinite one accepts every
    # disagreement. Both the value and the sigma are checked.
    @test _sigma_outcome(1.0, 2.0, 0.1, Inf; k=3).status === :error
    @test _sigma_outcome(1.0, 2.0, NaN, 0.1; k=3).status === :error
    @test _sigma_outcome(1.0, Inf, 0.1, 0.1; k=3).status === :error
    @test _sigma_outcome(1.0, 2.0, 0.1, 0.1; k=3).status === :fail
end

@testset "derivation_reach measures the other planes rather than guessing them" begin
    reach = derivation_reach()
    @test !isempty(reach)
    # A supplied temperature is not a second fetched value, so every reported hub
    # really does hold two atlas numbers the relation connects.  (`relations` is
    # non-empty by the producer's own filter, so asserting THAT would be
    # tautological; what is not free is that the hub implements enough.)
    @test all(r -> r.quantities ≥ 2, reach)
    @test any(r -> r.model === QAtlas.TFIM && r.bc === Infinite, reach)
end

# The hyperscaling suppression has no declared user today, so nothing would
# exercise it. It also cannot work route-by-route: with no dimension, `d` never
# enters the data, `consistency_report` never reaches Josephson, and there is no
# route to mark — the silent absence the exclusion exists to prevent. So the skip
# is emitted on its own, and both halves are asserted here.
@testset "a hub without a dimension says so, rather than losing the route quietly" begin
    free = QAtlas.SweptExponents(
        Universality{:Ising},
        Infinite,
        (d=[2],),
        nothing,
        Pair{Symbol,Vector{Symbol}}[],
        3.0,
        "an infinite-range hub has no spatial dimension",
        String[],
    )
    outcomes = Dict(
        c.id => QAtlas.run_generated_check(c) for c in QAtlas._scaling_checks(free)
    )
    for r in ("Josephson", "QuantumHyperscaling")
        id = "derivation/scaling/Universality{:Ising}/Infinite/" * r
        @test haskey(outcomes, id)
        @test outcomes[id].status === :skip
        @test occursin("hyperscaling does not apply", outcomes[id].detail)
    end
    # ...and the routes that need no dimension still judge, so dropping it costs
    # the hyperscaling laws and not the row.
    @test count(o -> o.status === :pass, values(outcomes)) >= 6
    # The same hub WITH a dimension reaches Josephson as a real route, which is
    # what makes the two ids above an exclusion rather than a restatement.
    withd = QAtlas.SweptExponents(
        Universality{:Ising},
        Infinite,
        (d=[2],),
        :sweep,
        Pair{Symbol,Vector{Symbol}}[],
        3.0,
        "",
        String[],
    )
    @test any(endswith("/ν/Josephson"), [c.id for c in QAtlas._scaling_checks(withd)])
end

@testset "k is the declaration's, not a constant in the generator" begin
    # Same table twice, differing only in k. The threshold is taken FROM the data
    # rather than guessed: run once at a k nothing can fail, read the largest
    # deviation the routes actually show, and straddle it.
    function outcomes(k)
        s = QAtlas.SweptExponents(
            Universality{:Ising},
            Infinite,
            (d=[3],),
            :sweep,
            Pair{Symbol,Vector{Symbol}}[],
            k,
            "",
            String[],
        )
        return [QAtlas.run_generated_check(c) for c in QAtlas._scaling_checks(s)]
    end
    wide = outcomes(1e3)
    @test !any(o -> o.status === :fail, wide)
    worst = maximum(o -> o.rel_err, wide)          # `rel_err` carries n_sigma
    @test worst > 0                                 # the table is not exact, so k can bite
    @test any(o -> o.status === :fail, outcomes(worst / 2))
    @test !any(o -> o.status === :fail, outcomes(worst * 2))
end
