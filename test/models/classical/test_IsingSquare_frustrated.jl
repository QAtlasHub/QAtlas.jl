# IsingSquare on a NON-BIPARTITE torus (#824): Z is not even in J there.
#
# The regression this guards was a `abs(K)` inside Kaufman's closed form, added
# with a correct-sounding justification — "on the bipartite square lattice the
# antiferromagnet maps to the ferromagnet, so Z is even in J".  True, and the
# periodic Lx × Ly lattice is bipartite only when BOTH sides are even.  With
# either side odd the wrap-around closes an odd cycle and the antiferromagnet is
# frustrated.
#
# It was invisible to every existing test because they all used even × even.  So
# the assertions here are chosen to need no stored reference and to be impossible
# to satisfy by accident:
#
#   * an EXHAUSTIVE enumeration over all 2^(Lx·Ly) configurations, which shares
#     nothing with either production route but the Boltzmann weight;
#   * the frustration INEQUALITY, which is a structural law, not a number;
#   * the torus symmetry Z(Lx, Ly) = Z(Ly, Lx), which exercises the axis swap.

using Test
using QAtlas
using QAtlas: IsingSquare, PartitionFunction, fetch, PBC
using QAtlas: FreeEnergy, Energy, ThermalEntropy, SpecificHeat

# Prefixed per file: test files share one `Main`, so an unprefixed `_Z` would be
# a global waiting to collide with whatever the shard planner puts next to it.
_isq_Z(Lx, Ly, J; β=1.0) = fetch(IsingSquare(), PartitionFunction(); β=β, Lx=Lx, Ly=Ly, J=J)

"""
    _isq_Z_bruteforce(Lx, Ly, β, J) -> Float64

Partition function of the `Lx × Ly` periodic square-lattice Ising model by summing
over all `2^(Lx·Ly)` spin configurations.

The point of doing it the stupid way: it shares no code path, no fermionisation
and no transfer matrix with either production route — only the definition of the
energy.  That is what makes it an oracle rather than a second opinion.
"""
function _isq_Z_bruteforce(Lx::Int, Ly::Int, β::Float64, J::Float64)
    n = Lx * Ly
    n <= 16 || error("_isq_Z_bruteforce: $n sites is too many")
    idx(x, y) = ((x % Lx) * Ly + (y % Ly)) + 1     # 0-based (x, y), PBC both ways
    Z = 0.0
    σ = Vector{Int}(undef, n)
    for m in 0:(2 ^ n - 1)
        @inbounds for i in 1:n
            σ[i] = ((m >> (i - 1)) & 1) == 1 ? 1 : -1
        end
        E = 0
        for x in 0:(Lx - 1), y in 0:(Ly - 1)
            s = σ[idx(x, y)]
            E += s * σ[idx(x + 1, y)]              # bond to the right
            E += s * σ[idx(x, y + 1)]              # bond upward
        end
        Z += exp(β * J * E)
    end
    return Z
end

@testset "brute force agrees on a FRUSTRATED torus, both signs of J" begin
    # 3×3 is the smallest non-bipartite torus; 2^9 = 512 configurations.
    for (Lx, Ly) in ((3, 3), (4, 3), (3, 4)), J in (-0.5, +0.5), β in (0.6, 1.0)
        Lx * Ly <= 16 || continue
        @test _isq_Z(Lx, Ly, J; β=β) ≈ _isq_Z_bruteforce(Lx, Ly, β, J) rtol = 1e-10
    end
    # …and on a bipartite one, where the old code was accidentally right
    for J in (-0.5, +0.5)
        @test _isq_Z(4, 4, J) ≈ _isq_Z_bruteforce(4, 4, 1.0, J) rtol = 1e-10
    end
    # the value #824 was filed on, pinned explicitly
    @test log(_isq_Z(3, 3, -0.5)) ≈ 7.8302303635 atol = 1e-9
    @test log(_isq_Z(3, 3, +0.5)) ≈ 9.9251503709 atol = 1e-9
end

@testset "frustration is an INEQUALITY, so it needs no reference data" begin
    # An odd cycle cannot be satisfied by any assignment, so the antiferromagnet
    # carries strictly less weight than the ferromagnet. This assertion cannot
    # rot: no stored number, no tolerance on a value, just a strict order.
    for (Lx, Ly) in ((3, 3), (4, 3), (3, 4), (5, 3), (3, 5), (5, 5), (3, 8), (7, 4))
        @test _isq_Z(Lx, Ly, -0.5) < _isq_Z(Lx, Ly, +0.5)
    end
    # and it must NOT hold where the lattice is bipartite — equality, exactly,
    # since the two are related by a gauge transformation rather than a limit
    for (Lx, Ly) in ((4, 4), (6, 4), (4, 6), (8, 8), (2, 6))
        @test _isq_Z(Lx, Ly, -0.5) == _isq_Z(Lx, Ly, +0.5)
    end
    # colder ⇒ more frustration ⇒ a wider gap
    gap(β) = log(_isq_Z(3, 3, +0.5; β=β)) - log(_isq_Z(3, 3, -0.5; β=β))
    @test gap(0.3) < gap(0.6) < gap(1.0) < gap(2.0)
end

@testset "the DERIVATIVE quantities are right too, not merely non-throwing" begin
    # `_ising_sq_log_z` feeds FreeEnergy, Energy, ThermalEntropy and SpecificHeat,
    # so the #824 sign loss was never confined to PartitionFunction. The frustrated
    # branch also has to keep propagating ForwardDiff Duals, which its docstring
    # requires and which the transfer-matrix route is not obviously able to do —
    # so this checks the VALUES against finite differences of the brute-force
    # log Z, not just that nothing throws.
    logZbf(Lx, Ly, β, J) = log(_isq_Z_bruteforce(Lx, Ly, β, J))
    for (Lx, Ly, J) in ((3, 3, -0.5), (3, 3, +0.5), (4, 3, -0.5))
        N, β, h = Lx * Ly, 1.0, 1e-5
        f(b) = logZbf(Lx, Ly, b, J)
        e_bf = -(f(β + h) - f(β - h)) / (2h) / N
        c_bf = β^2 * (f(β + h) - 2f(β) + f(β - h)) / h^2 / N
        e_at = fetch(IsingSquare(), Energy(), PBC(); beta=β, Lx=Lx, Ly=Ly, J=J)
        c_at = fetch(IsingSquare(), SpecificHeat(), PBC(); beta=β, Lx=Lx, Ly=Ly, J=J)
        @test e_at ≈ e_bf atol = 1e-7
        @test c_at ≈ c_bf atol = 1e-4        # finite-difference noise dominates here
    end
    # free energy is just -log Z / (β N), so it inherits the fix directly
    @test fetch(IsingSquare(), FreeEnergy(), PBC(); beta=1.0, Lx=3, Ly=3, J=-0.5) ≈
        -log(_isq_Z(3, 3, -0.5)) / 9 rtol = 1e-12
end

@testset "the torus is symmetric under swapping its axes" begin
    # Z(Lx, Ly) = Z(Ly, Lx). This is what exercises the axis swap in the
    # frustrated route, which puts the exponential cost on the SHORTER side —
    # so 3×14 must work and must equal 14×3.
    for (Lx, Ly) in ((3, 4), (5, 3), (3, 14), (7, 4)), J in (-0.5, +0.5)
        @test _isq_Z(Lx, Ly, J) ≈ _isq_Z(Ly, Lx, J) rtol = 1e-12
    end
end

@testset "it refuses rather than silently returning the ferromagnetic value" begin
    # Beyond the transfer-matrix cap on the SHORT side there is no exact route,
    # and the closed form does not apply. Refusing is the point of #824: the bug
    # was not a missing feature, it was a wrong number where a refusal belonged.
    # NON-bipartite with the short side over the cap.  Note 14×14 does NOT throw
    # and must not: both sides even, so it folds onto |K| legitimately — which is
    # exactly the trap this testset nearly fell into.
    @test_throws ArgumentError _isq_Z(13, 13, -0.5)
    @test_throws ArgumentError _isq_Z(15, 14, -0.5)
    err = try
        _isq_Z(13, 13, -0.5)
    catch e
        e
    end
    @test occursin("FRUSTRATED", err.msg)
    # the same sizes are fine ferromagnetically …
    @test isfinite(_isq_Z(13, 13, +0.5))
    @test isfinite(_isq_Z(15, 14, +0.5))
    # … and a BIPARTITE torus of the same scale answers for either sign, because
    # there the fold is a gauge transformation rather than an approximation
    @test _isq_Z(14, 14, -0.5) == _isq_Z(14, 14, +0.5)
    @test isfinite(_isq_Z(14, 14, -0.5))
end

@testset "the closed form now states its domain instead of folding silently" begin
    # `_ising_sq_log_z_torus` is Kaufman's formula and is written for K ≥ 0.
    # It used to accept any K and quietly use |K|; now it says so.
    @test isfinite(QAtlas._ising_sq_log_z_torus(3, 3, 0.5))
    @test_throws ArgumentError QAtlas._ising_sq_log_z_torus(3, 3, -0.5)
    @test QAtlas._ising_sq_bipartite_torus(4, 4)
    @test !QAtlas._ising_sq_bipartite_torus(3, 4)
    @test !QAtlas._ising_sq_bipartite_torus(4, 3)
    @test !QAtlas._ising_sq_bipartite_torus(3, 3)
end
