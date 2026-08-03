# The `:quench_relaxes_to_gge` identity edge — lim_{t→∞} ⟨σˣ⟩(t) = ⟨σˣ⟩_GGE.
#
# The generated checks themselves run in test/generated/test_identity_rest.jl,
# which picks up any new edge automatically.  What needs its own file is the set
# of CLAIMS the edge's header makes: that it discriminates, that `t = 800` is a
# plateau rather than a step on the way to something better, and that its
# independence has a stated LIMIT.

using Test
using QAtlas
using QAtlas: TFIM, Infinite, fetch, QuenchLocalMagnetization, GGEValue, MagnetizationX
using QAtlas: RenyiEntropy, VonNeumannEntropy, _quantity_instance, generated_checks

_Q(mf, m0, t) = fetch(mf, QuenchLocalMagnetization{:x}(), Infinite(); initial=m0, t=t)
_G(mf, m0) = fetch(mf, GGEValue(MagnetizationX()), Infinite(); initial=m0)

@testset "the edge materialises, and only on the hub that has both" begin
    cs = filter(
        c -> occursin("quench_relaxes_to_gge", c.id), generated_checks(; kinds=(:identity,))
    )
    @test length(cs) == 4                       # 4 initial fields × 1 time
    @test all(c -> occursin("/TFIM/", c.id), cs)
    # Infinite only: `GGEValue{MagnetizationX}` is registered nowhere else, and
    # the OBC quench route needs a site index the GGE has no analogue for.
    @test all(c -> occursin("/Infinite/", c.id), cs)
    for c in cs
        r = QAtlas.run_generated_check(c)
        @test r.status === :pass
        @test r.abs_err < 1e-6                  # measured worst case 2.2e-7
    end
end

@testset "it DISCRIMINATES — the tolerance is not just larger than the residual" begin
    # A check that passes is worth nothing unless the wrong answer fails it.
    # Compare the quench from h₀ against the GGE for a DIFFERENT h₀, at the same
    # final field, and require the gap to dwarf the tolerance.
    mf = TFIM(1.0, 1.0)
    tol = 1e-5
    for h0 in (0.2, 0.5, 0.8, 1.5, 2.0, 3.0)
        m0 = TFIM(1.0, h0)
        resid = abs(_Q(mf, m0, 800.0) - _G(mf, m0))
        margin = abs(_Q(mf, m0, 800.0) - _G(mf, TFIM(1.0, h0 + 0.1)))
        @test resid < tol / 10                  # ≥10× headroom below the tolerance
        @test margin > tol * 100                # ≥100× above it
        @test margin / max(resid, eps()) > 1e3  # and ≥10³ separation between them
    end
end

@testset "t = 800 is a PLATEAU, not a step toward something better" begin
    # The header claims raising `t` does not help, because the cos(2Λt) term
    # eventually oscillates faster than the quadrature can follow.  If that were
    # wrong — if the residual kept falling — the honest edge would use a larger
    # `t` and a tighter tolerance, so the claim is load-bearing and pinned here.
    mf, m0 = TFIM(1.0, 1.0), TFIM(1.0, 0.5)
    g = _G(mf, m0)
    r(t) = abs(_Q(mf, m0, t) - g)
    r100, r800, r3200 = r(100.0), r(800.0), r(3200.0)
    @test r100 > 10 * r800                      # it DOES improve up to ~800 …
    @test r3200 > r800 / 10                     # … and then stops: same order, no gain
    @test r3200 < 10 * r800
end

@testset "the independence has a stated LIMIT: a shared convention" begin
    # The two integrands are different expressions, which is what makes the edge
    # a real check.  But they reach the Bogoliubov angle through SEPARATE helper
    # families that are duplicate implementations of one convention — so the edge
    # is blind to an error in that convention, both sides moving together.
    #
    # Asserting the duplication is the only way that caveat stays true: if the
    # helpers ever diverge, this fails and the header's scope statement must be
    # rewritten rather than silently becoming wrong.
    for k in (0.3, 1.0, 2.5), h in (0.5, 1.5, 2.0)
        @test QAtlas._tfim_two_theta(h, 1.0, k) ≈ QAtlas._tfim_gge_two_theta(k, 1.0, h) rtol =
            1e-12
        @test QAtlas._tfim_lambda(h, 1.0, k) ≈ QAtlas._tfim_dispersion(k, 1.0, h) rtol =
            1e-12
    end
end

@testset "_quantity_instance: the TYPE determines the instance, or it throws" begin
    # `GGEValue{Q}` is a functor over quantities (#819 item 3) and carries a
    # field, so the old rule ("zero fields or error") made it unreferenceable
    # from any constraint edge.  The new rule is narrower than "wrappers are OK".
    @test _quantity_instance(GGEValue{MagnetizationX}) == GGEValue(MagnetizationX())
    @test _quantity_instance(VonNeumannEntropy) == VonNeumannEntropy()

    # …and a field the type does NOT determine still throws: no type names which
    # α is meant, so an edge over `RenyiEntropy` has to say.
    @test_throws ArgumentError _quantity_instance(RenyiEntropy)
    err = try
        _quantity_instance(RenyiEntropy)
    catch e
        e
    end
    @test occursin("does not determine", err.msg)
end
