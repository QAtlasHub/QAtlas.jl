# ─────────────────────────────────────────────────────────────────────────────
# Verification: nearest-neighbor tight-binding on the honeycomb lattice
#
# Cross-validate QAtlas's closed-form Bloch spectrum against direct
# diagonalization of the real-space tight-binding Hamiltonian built from
# `Lattice2D.bonds(lat)`.
#
# Test sizes: 2×2 through 4×4 honeycomb PBC.
#
# Additional structural checks motivated by the chiral (sublattice)
# symmetry of bipartite honeycomb:
#   * Σ λᵢ = 0   (tr H = 0)
#   * {λᵢ} is symmetric about 0  (λ ↔ −λ)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

const T_HOP = 1.0

@testset "Graphene TB — real-space vs Bloch closed form" begin
    for (Lx, Ly) in [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]
        lat = build_lattice(Honeycomb, Lx, Ly)

        @testset "$(Lx)×$(Ly) honeycomb PBC" begin
            H = build_tight_binding(lat, T_HOP)
            λ_real = sort(eigvals(Symmetric(H)))
            λ_bloch = QAtlas.fetch(
                Graphene(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=T_HOP
            )

            @test length(λ_real) == 2 * Lx * Ly
            @test length(λ_bloch) == 2 * Lx * Ly
            @test λ_real ≈ λ_bloch atol = 1e-10

            # Chiral (sublattice) symmetry: spectrum symmetric about 0
            @test sum(λ_real) ≈ 0 atol = 1e-10
            @test sort(λ_real) ≈ sort(-λ_real) atol = 1e-10
        end
    end

    @testset "2×2 honeycomb — hand-computable spectrum" begin
        # k_{mn} = (2π m/2, 2π n/2) = (πm, πn) → spectrum {-3, -1, -1, -1, 1, 1, 1, 3}
        λ = QAtlas.fetch(Graphene(), TightBindingSpectrum(); Lx=2, Ly=2, t=1.0)
        @test λ ≈ [-3.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 3.0] atol = 1e-12
    end

    @testset "3×3 honeycomb — Dirac points give zero modes" begin
        # At (m,n) = (1,2) and (2,1) (the K, K' points of the finite BZ),
        # |f(k)|² = 0 exactly → four zero modes (two k-points × two bands).
        λ = QAtlas.fetch(Graphene(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
        zero_modes = count(x -> abs(x) < 1e-10, λ)
        @test zero_modes == 4
    end

    @testset "t scaling: λ(t) = t · λ(1)" begin
        λ1 = QAtlas.fetch(Graphene(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
        for t in (0.5, 2.0, 3.7)
            λt = QAtlas.fetch(Graphene(), TightBindingSpectrum(); Lx=3, Ly=3, t=t)
            @test λt ≈ t .* λ1 rtol = 1e-12
        end
    end
end
