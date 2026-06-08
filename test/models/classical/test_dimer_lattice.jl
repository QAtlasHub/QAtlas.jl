# ─────────────────────────────────────────────────────────────────────────────
# Test: DimerLattice — close-packed dimers on the square lattice (Kasteleyn).
#
# The rigorous net is the INDEPENDENT brute-force perfect-matching enumeration:
# the Kasteleyn-Temperley-Fisher closed form must equal a direct recursive count
# for every small grid (an unarguable witness, validated itself against known
# literature counts).  The thermodynamic entropy G/π is then checked by the
# finite-grid (ln Z)/N limit.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using QAtlas:
    DimerLattice, PartitionFunction, ResidualEntropy, FreeEnergy, OBC, Infinite, fetch
using Base.MathConstants: catalan

# Independent witness: brute-force perfect-matching count by recursive
# backtracking (no reference to the KTF formula).
function _brute_dimer_count(m::Int, n::Int)
    iseven(m * n) || return 0
    occ = falses(m, n)
    function rec()
        i0 = 0
        j0 = 0
        @inbounds for i in 1:m, j in 1:n
            if !occ[i, j]
                i0, j0 = i, j
                break
            end
        end
        i0 == 0 && return 1
        t = 0
        if j0 < n && !occ[i0, j0 + 1]            # place horizontally
            occ[i0, j0] = occ[i0, j0 + 1] = true
            t += rec()
            occ[i0, j0] = occ[i0, j0 + 1] = false
        end
        if i0 < m && !occ[i0 + 1, j0]            # place vertically
            occ[i0, j0] = occ[i0 + 1, j0] = true
            t += rec()
            occ[i0, j0] = occ[i0 + 1, j0] = false
        end
        return t
    end
    return rec()
end

@testset "DimerLattice — structural / error guards" begin
    # pin the brute-force witness itself against known literature counts
    @test _brute_dimer_count(2, 2) == 2
    @test _brute_dimer_count(2, 4) == 5            # 2×n is Fibonacci(n+1)
    @test _brute_dimer_count(4, 4) == 36
    @test _brute_dimer_count(6, 6) == 6728
    @test _brute_dimer_count(1, 2) == 1            # 1×n strip: 1 tiling for even n
    @test _brute_dimer_count(1, 4) == 1
    @test _brute_dimer_count(3, 4) == 11           # odd×even shape anchor

    # odd number of sites ⇒ no perfect matching
    @test fetch(DimerLattice(), PartitionFunction(); Lx=3, Ly=3) == 0
    @test fetch(DimerLattice(), PartitionFunction(); Lx=5, Ly=3) == 0

    # missing size errors loudly; struct vs kwargs agree
    @test_throws ErrorException fetch(DimerLattice(), PartitionFunction())
    @test fetch(DimerLattice(; Lx=4, Ly=4), PartitionFunction()) ==
        fetch(DimerLattice(), PartitionFunction(); Lx=4, Ly=4)

    # negative sizes rejected at construction; counts beyond exact Float64 integer
    # range (2^53) error loudly rather than return a silent Inf / wrong integer
    @test_throws ArgumentError DimerLattice(; Lx=-2, Ly=3)
    @test_throws ErrorException fetch(DimerLattice(), PartitionFunction(); Lx=14, Ly=14)

    # FreeEnergy density = −ResidualEntropy (unit weights ⇒ zero internal energy)
    @test fetch(DimerLattice(), FreeEnergy(), Infinite()) ==
        -fetch(DimerLattice(), ResidualEntropy(), Infinite())
end

# ── INDEPENDENT: KTF closed form == brute-force enumeration over a grid sweep ──
@testset "DimerLattice — KTF count == brute-force enumeration (m,n sweep)" begin
    for m in 1:6, n in 1:6
        @test fetch(DimerLattice(), PartitionFunction(); Lx=m, Ly=n) ==
            _brute_dimer_count(m, n)
    end
end

# ── ResidualEntropy: finite-grid (1/N) ln Z → G/π (Catalan/π) ─────────────────
@testset "DimerLattice — entropy per site → G/π (Catalan/π)" begin
    s∞ = catalan / π
    @test fetch(DimerLattice(), ResidualEntropy(), Infinite()) == s∞
    # (ln Z)/N from the KTF log-sum directly (avoids the count overflowing Float64
    # at large L); the count itself is validated against brute force above.
    sden(L) =
        sum(
            log(4 * cos(j * π / (L + 1))^2 + 4 * cos(k * π / (L + 1))^2) for
            j in 1:L, k in 1:L
        ) / (4 * L * L)
    @test abs(sden(40) - s∞) < abs(sden(10) - s∞)      # monotone approach
    @test abs(sden(80) - s∞) < abs(sden(40) - s∞)
    @test isapprox(sden(120), s∞; atol=5e-3)           # close at large L
end

# ── Verification cards (WHY-correct plane) ────────────────────────────────────
@testset "DimerLattice — verification cards" begin
    for (m, n) in ((2, 2), (4, 4), (6, 6), (2, 6), (4, 6))
        verify(
            DimerLattice(; Lx=m, Ly=n),
            PartitionFunction(),
            OBC();
            route=:ed_finite_size,
            fetch_kw=(; Lx=m, Ly=n),
            independent=Float64(_brute_dimer_count(m, n)),
            agree_within=1e-6,
            refs=[
                "Kasteleyn 1961: KTF Pfaffian count vs brute-force perfect-matching enumeration",
            ],
        )
    end
end
