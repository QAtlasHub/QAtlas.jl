# ─────────────────────────────────────────────────────────────────────────────
# Verification: NN tight-binding on the Union Jack (centred square) lattice
#
# Two sublattices: A (corner, 8-coordinated) and B (body-centre,
# 4-coordinated). A–A square-lattice bonds exist → NOT bipartite → no
# chiral symmetry and no flat-band guarantee.
#
# Zero src code — generic bloch_tb_spectrum only.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

const T_HOP = 1.0

@testset "UnionJack TB — generic Bloch vs real-space ED" begin
    for (Lx, Ly) in [(3, 3), (3, 4), (4, 4)]
        lat = build_lattice(UnionJack, Lx, Ly)

        @testset "$(Lx)×$(Ly) UnionJack PBC" begin
            H = build_tight_binding(lat, T_HOP)
            λ_real = sort(eigvals(Symmetric(H)))
            λ_bloch = bloch_tb_spectrum(UnionJack, Lx, Ly, T_HOP)

            @test length(λ_real) == 2 * Lx * Ly
            @test length(λ_bloch) == 2 * Lx * Ly
            @test λ_real ≈ λ_bloch atol = 1e-10

            # Structural: tr H = 0 (NN TB)
            @test tr(H) ≈ 0 atol = 1e-10

            # NOT bipartite: spectrum is NOT symmetric about zero
            @test !(sort(λ_real) ≈ sort(-λ_real))
        end
    end

    @testset "t scaling" begin
        λ1 = bloch_tb_spectrum(UnionJack, 3, 3, 1.0)
        for t in (0.5, 2.0)
            λt = bloch_tb_spectrum(UnionJack, 3, 3, t)
            @test λt ≈ t .* λ1 rtol = 1e-12
        end
    end
end
