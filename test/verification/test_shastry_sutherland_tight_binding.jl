# ─────────────────────────────────────────────────────────────────────────────
# Verification: NN tight-binding on the Shastry–Sutherland lattice
#
# Four sublattices per unit cell (square + dimer bonds). The unit cell
# has two bond types (type 1 = NN square, type 2 = dimer diagonal), but
# the generic builder treats all bonds with uniform hopping -t.
#
# The Shastry–Sutherland lattice is famous for its orthogonal-dimer
# structure; its tight-binding spectrum reflects the interplay between
# square-lattice connectivity and extra dimer hoppings.
#
# Zero src code — generic bloch_tb_spectrum only.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

const T_HOP = 1.0

@testset "Lattice2D.ShastrySutherland TB — generic Bloch vs real-space ED" begin
    for (Lx, Ly) in [(2, 2), (3, 3), (3, 4), (4, 4)]
        lat = build_lattice(Lattice2D.ShastrySutherland, Lx, Ly)

        @testset "$(Lx)×$(Ly) SS PBC" begin
            H = build_tight_binding(lat, T_HOP)
            λ_real = sort(eigvals(Symmetric(H)))
            λ_bloch = bloch_tb_spectrum(Lattice2D.ShastrySutherland, Lx, Ly, T_HOP)

            @test length(λ_real) == 4 * Lx * Ly
            @test length(λ_bloch) == 4 * Lx * Ly
            @test λ_real ≈ λ_bloch atol = 1e-10

            # Structural: tr H = 0 (NN TB, no on-site)
            @test tr(H) ≈ 0 atol = 1e-10
        end
    end

    @testset "t scaling" begin
        λ1 = bloch_tb_spectrum(Lattice2D.ShastrySutherland, 3, 3, 1.0)
        for t in (0.5, 2.0)
            λt = bloch_tb_spectrum(Lattice2D.ShastrySutherland, 3, 3, t)
            @test λt ≈ t .* λ1 rtol = 1e-12
        end
    end
end
