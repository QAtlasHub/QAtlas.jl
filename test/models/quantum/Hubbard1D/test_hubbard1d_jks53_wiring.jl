# Hubbard1D/FreeEnergy/Infinite routed to the eq (53) solver (#798).
#
# WHAT CHANGED.  The row used to call the eq (47) NLIE, which carried four
# defects at once: it convolved the c channel with `K1` where eq (47) says
# `K1bar`, used `log(1+b)` where the paper says `log(1+1/b)`, kept one boundary
# value per function where eq (51) needs two, and — the one that dominated the
# error — passed `μ = U/2` straight through as though it were the PAPER's half
# filling.  JKS write the interaction symmetrically, `U(n↓−½)(n↑−½)`, so their
# half filling is `μ = 0`.  The old route was solving a DOPED system and
# comparing it against half-filled ED.
#
# THIS IS AN IMPROVEMENT, NOT A CLOSURE.  eq (53) is better everywhere and
# converges where the old route returned `NaN`, but it still drifts at low
# temperature, and that drift is NOT a resolution artifact — measured, it is
# converged in the solver's own knobs:
#
#     β = 1, Nw/Nn/x_max = 96/48/32   →  −2.669123
#                          160/80/32  →  −2.668892
#                          256/128/64 →  −2.668871
#
# while ED gives −3.125427.  So the remaining gap is in the formulation or its
# implementation, and #798 stays open.  The assertions below therefore BOUND the
# known error from both sides rather than asserting agreement: a silent
# regression fails, and so does a silent fix — the latter on purpose, because a
# fix that nobody notices is a fix nobody documents.

using Test
using LinearAlgebra, SparseArrays
using QAtlas
using QAtlas: Hubbard1D, FreeEnergy, Infinite, fetch

include(joinpath(@__DIR__, "..", "..", "..", "util", "hubbard_ed.jl"))

const _JKS53_U = 4.0
const _JKS53_T = 1.0
const _JKS53_MU = 2.0            # plain-convention half filling, μ = U/2

# Memoised: each NLIE solve is ~10 s and the same β is asserted on from several
# testsets. Caching the RESULT changes no assertion — it removes repeated solves
# of an identical problem, taking the file from 5m30s to well under a minute.
const _JKS53_CACHE = Dict{Any,Float64}()
const _JKS53_ED_CACHE = Dict{Any,Float64}()

function _jks53_f(beta; kw...)
    key = (beta, sort!(collect(pairs(kw)); by=first))
    return get!(_JKS53_CACHE, key) do
        return fetch(
            Hubbard1D(; t=_JKS53_T, U=_JKS53_U, μ=_JKS53_MU),
            FreeEnergy(),
            Infinite();
            beta=beta,
            kw...,
        )
    end
end
function _jks53_ed(beta; N=6)
    return get!(_JKS53_ED_CACHE, (beta, N)) do
        return _ed_hubbard_free_energy(N, _JKS53_T, _JKS53_U, _JKS53_MU, beta; pbc=true)
    end
end

@testset "the convention map is executable, not prose" begin
    # f_paper(μ) = f_plain(μ + U/2) + U/4  ⇒  f_plain(μ) = f_paper(μ − U/2) − U/4
    U = _JKS53_U
    @test QAtlas._jks_paper_mu(U / 2, U) == 0.0        # plain half filling ↦ paper's
    @test QAtlas._jks_paper_mu(0.0, U) == -U / 2
    @test QAtlas._jks_plain_offset(U) == U / 4
    # the round trip: going to the paper convention and back is the identity on μ
    for mu in (0.0, 1.0, 2.0, 3.5)
        @test QAtlas._jks_paper_mu(mu, U) + U / 2 ≈ mu
    end
    # …and it is the mapping the wired fetch actually uses, not a parallel copy
    v_direct =
        QAtlas.Hubbard1DJKSNLIE.hubbard1d_jks53_free_energy(
            _JKS53_T, _JKS53_U, QAtlas._jks_paper_mu(_JKS53_MU, U), 0.1
        ) - QAtlas._jks_plain_offset(U)
    @test _jks53_f(0.1) ≈ v_direct rtol = 1e-12
end

@testset "the high-T limit is reproduced, and it needs no reference data" begin
    # Four states per site, all equally weighted as β → 0, so f → −ln4/β with no
    # free parameters at all — the cheapest possible independent check, and the
    # one the old route failed by ≈ 3 at β = 0.1.
    for b in (0.01, 0.02, 0.05)
        @test _jks53_f(b) ≈ -log(4) / b rtol = 0.04
    end
    # and the approach is from below in |f| — the interaction and hopping
    # corrections are O(1) against a term diverging as 1/β
    @test abs(_jks53_f(0.01) + log(4) / 0.01) < abs(_jks53_f(0.05) + log(4) / 0.05)
end

@testset "agreement with ED, and the ED oracle is itself converged" begin
    # Assert the oracle before leaning on it: at these β the finite-N free
    # energy has stopped moving, so a disagreement is the NLIE's, not the box's.
    for b in (0.1, 0.3, 0.5)
        @test _jks53_ed(b; N=5) ≈ _jks53_ed(b; N=6) atol = 1e-3
    end
    # high T: the regime the row is registered for
    for b in (0.05, 0.1)
        @test _jks53_f(b) ≈ _jks53_ed(b) rtol = 0.01
    end
end

@testset "β ≥ 0.3 returns a number instead of NaN" begin
    # The concrete #798 win: the old route's β-continuation collapsed and `fetch`
    # returned NaN from β ≈ 0.3 upward.
    for b in (0.3, 0.5, 1.0)
        v = _jks53_f(b)
        @test isfinite(v)
        @test v < 0                                  # a free energy, not a sign error
    end
    # monotone in β, as a free energy of a stable system must be
    vs = [_jks53_f(b) for b in (0.1, 0.2, 0.3, 0.5, 1.0)]
    @test issorted(vs)                               # −14.9 < −8.0 < … < −2.7
end

@testset "the REMAINING low-T error is bounded, in both directions" begin
    # Not "it agrees". It does not, and pretending otherwise is how a known
    # defect becomes an unknown one. Lower bounds fail if someone fixes this
    # without updating the docs; upper bounds fail if it regresses.
    err(b) = abs(_jks53_f(b) - _jks53_ed(b)) / abs(_jks53_ed(b))
    @test 0.0005 < err(0.1) < 0.02        # measured 0.52 %
    @test 0.01 < err(0.3) < 0.06          # measured 3.7 %
    @test 0.08 < err(1.0) < 0.25          # measured 14.6 %
    # it grows with β — which is the shape that says "low-T formulation", not
    # "noise", and is why the row carries a valid_domain rather than a tolerance
    @test err(0.1) < err(0.3) < err(1.0)
end

@testset "the residual error is NOT the discretisation" begin
    # If refining the solver's own grids moved the answer toward ED, this would
    # be a knob problem and the row would just need better defaults. It does not:
    # the answer is converged to ~2e-4 while sitting 0.46 away from ED.
    coarse = _jks53_f(1.0; Nw=96, Nn=48, x_max=32.0)
    fine = _jks53_f(1.0; Nw=160, Nn=80, x_max=32.0)
    @test coarse ≈ fine atol = 1e-3
    @test abs(coarse - _jks53_ed(1.0)) > 0.3
end
