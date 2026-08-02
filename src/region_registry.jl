# region_registry.jl — declared REGION edges (core/region_checks.jl).
#
# The entanglement inequalities AbstractQAtlas carries but that no generator on
# this side could reach (#780).  Their slots are the same quantity on different
# REGIONS, not different quantity types, so `@bound_edge` rightly rejects them.

# ── the entropy inequalities on three adjacent blocks ─────────────────
#
# WHY ADJACENT, and why 2-site blocks.  Every union `region_report` forms from
# `A = 1:2`, `B = 3:4`, `C = 5:6` — `A∪B = 1:4`, `B∪C = 3:6`, `A∪B∪C = 1:6` — is
# a single contiguous interval.  That is the only shape a free-fermion hub can
# answer for: outside one interval the Jordan-Wigner string does not factorise
# and the covariance submatrix returns the FERMIONIC entropy, which differs from
# the spin one by 0.5-0.6 nats and which no entropy inequality would flag (#783).
# Dense-ED hubs take any region (#786), but sharing the blocks keeps every hub on
# the same instances, so a disagreement between two hubs means physics and not a
# difference of setup.
#
# `A∪C = {1,2,5,6}` is the family's ONE non-contiguous member (the family is
# every union of blocks).  Dense-ED hubs answer it and gain the instances it
# supports; free-fermion hubs refuse it, so their bag has 6 regions instead of 7.
# MEASURED: TFIM/OBC gives |bag| = 6, XXZ1D/OBC and S1Heisenberg1D/OBC give 7 --
# the Jordan-Wigner guard from #783 doing exactly its job.
#
# WHICH instances that costs is not obvious, and two plausible readings of it are
# both wrong; the test asserts the right one in BOTH directions.  Of TFIM/OBC's
# 18 instances, 8 skip, for TWO different reasons:
#
#   * the region ITSELF is non-contiguous -- the bipartite pair
#     (A = {1,2,5,6}, B = {3,4}) needs S(A), and is refused for A even though
#     A ∪ B = {1..6} is contiguous;
#   * a UNION is non-contiguous -- the triples put B at an end, e.g.
#     (A = {1,2}, B = {5,6}, C = {3,4}) needs A ∪ B = {1,2,5,6}.
#
# Note what is NOT a reason: none of these relations forms A ∪ C.  The bipartite
# pair needs A ∪ B; both triples are stated over A ∪ B and B ∪ C (strong
# subadditivity adds A ∪ B ∪ C and B, weak monotonicity adds A and C).  So the
# in-order triple (A = 1:2, B = 3:4, C = 5:6) needs only contiguous unions and
# does NOT skip, even though it names both {1,2} and {5,6}.
#
# N = 8 leaves sites 7-8 outside every region, so `A∪B∪C` is a proper subsystem
# and its entropy is not the trivial 0 of a pure whole.
# MEASURED COST, and the one exclusion it forces.  One bag (7 region entropies)
# at N = 8, beta = Inf:
#
#     TFIM/OBC              0.0 s      XXZ1D/OBC          0.38 s
#     Heisenberg1D/OBC      0.67 s     S1Heisenberg1D/OBC 338.7 s
#
# S1Heisenberg1D is three orders of magnitude off the rest.  It is not the
# spin-1 dimension alone (3^8 = 6561 against 2^8 = 256): its `fetch` builds the
# FULL 6561x6561 thermal rho even at beta = Inf, where the state is pure and the
# reduced state could come from the 6561-element VECTOR instead.  That is a fixable
# implementation cost, not a property of the physics, so it is excluded here with
# the measurement recorded rather than worked around by shrinking the blocks for
# everyone -- see the follow-up issue.
const _REGION_COST_EXCLUSIONS = [
    S1Heisenberg1D =>
        "one region bag costs 338.7 s at N = 8 (MEASURED), against " *
        "0.0-0.7 s for every other hub: its fetch builds the full " *
        "3^N x 3^N thermal rho even at beta = Inf, where the state is " *
        "pure; re-enable once that path uses the state vector",
]

@region(
    :entanglement_regions,
    blocks = [[1, 2], [3, 4], [5, 6]],
    finite_N = 8,
    sweep = (beta=[Inf],),
    exclusions = _REGION_COST_EXCLUSIONS,
    # Entropy inequalities are exact statements about a state, not numerical
    # agreements between two routes, so the tolerance is a round-off guard and
    # nothing more.  A violation is a broken reduced state, at any magnitude.
    atol = 1e-9,
    notes =
        "Subadditivity, Araki-Lieb, weak monotonicity and strong " *
        "subadditivity over three adjacent blocks, auto-discovered by " *
        "AbstractQAtlas's region_report (#780).",
)
