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

@testset "exponent_sweep! refuses the declarations that would read as coverage" begin
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Ising}, Infinite; sweep=(d=[2],)
    )                                                  # duplicate (model, bc, sweep)
    @test_throws ArgumentError exponent_sweep!(Universality{:Potts3}, PBC; k=0)
    @test_throws ArgumentError exponent_sweep!(Universality{:Potts3}, PBC; sweep=(;))
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Potts3}, PBC; dimension=:nope, sweep=(d=[2],)
    )
    # dimension=nothing suppresses the hyperscaling routes; unexplained, that is
    # a silent exclusion wearing the shape of a pass.
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Potts3}, PBC; dimension=nothing
    )
    @test_throws ArgumentError exponent_sweep!(
        Universality{:Potts3}, PBC; refused="because", notes="also because"
    )
    # ...and nothing above left a row behind.
    @test !any(s -> s.bc === PBC, EXPONENT_SWEEPS)
end

@testset "derivation_reach measures the other planes rather than guessing them" begin
    reach = derivation_reach()
    @test !isempty(reach)
    # A supplied temperature is not a second fetched value, so every reported hub
    # really does hold two atlas numbers the relation connects.
    @test all(r -> !isempty(r.relations), reach)
    @test all(r -> r.quantities ≥ 2, reach)
    @test any(r -> r.model === QAtlas.TFIM && r.bc === Infinite, reach)
end
