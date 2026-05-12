# ─────────────────────────────────────────────────────────────────────────────
# Verification: nearest-neighbor tight-binding on the triangular lattice
#
# Cross-validate QAtlas's closed-form single-band Bloch spectrum against
# direct diagonalization of the real-space tight-binding Hamiltonian
# built from `Lattice2D.bonds(lat)`.
#
# Test sizes: 3×3 through 6×6 triangular PBC. Sizes with Lx, Ly multiple
# of 3 have the K-points in the discrete BZ, reproducing the hallmark
# asymmetric band edge -6t / +3t.
#
# Additional structural checks:
#   * tr H = 0                          (NN TB has no on-site energy)
#   * Band edges -6t at Γ (unique) and +3t at K / K' (degeneracy 2)
#     when both Lx and Ly are divisible by 3
#   * Spectrum is NOT symmetric about zero (frustration, non-bipartite)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

const T_HOP = 1.0

@testset "Triangular TB — real-space vs Bloch closed form" begin
    for (Lx, Ly) in [(3, 3), (3, 4), (4, 4), (5, 5), (6, 6)]
        lat = build_lattice(Triangular, Lx, Ly)

        @testset "$(Lx)×$(Ly) triangular PBC" begin
            H = build_tight_binding(lat, T_HOP)
            λ_real = sort(eigvals(Symmetric(H)))
            λ_bloch = QAtlas.fetch(
                QAtlas.Triangular(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=T_HOP
            )

            @test length(λ_real) == Lx * Ly
            @test length(λ_bloch) == Lx * Ly
            @test λ_real ≈ λ_bloch atol = 1e-10

            # Structural: tr H = 0
            @test tr(H) ≈ 0 atol = 1e-10

            # Γ-point gives -6t unique minimum
            @test λ_real[1] ≈ -6.0 atol = 1e-10
            @test count(x -> abs(x - (-6.0)) < 1e-10, λ_real) == 1
        end
    end

    @testset "3×3 hand-computable spectrum" begin
        # Allowed (m, n) ∈ {0, 1, 2}²:
        #   Γ = (0,0)                              → -6
        #   (1,1), (2,2)                           →  0 (θ₂−θ₁ = 0)
        #   (1,0), (0,1), (2,0), (0,2)             →  0 (one θ = 0)
        #   (1,2)                                  → +3 (K)
        #   (2,1)                                  → +3 (K′)
        # Sorted: [-6, 0×6, 3, 3]
        λ = QAtlas.fetch(QAtlas.Triangular(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
        @test λ ≈ [-6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 3.0] atol = 1e-12
    end

    @testset "K / K' degeneracy at multiples of 3" begin
        # When both Lx and Ly are divisible by 3 the K-points
        # (m/Lx, n/Ly) = (1/3, 2/3) and (2/3, 1/3) are allowed momenta
        # and contribute +3t eigenvalues.
        for (Lx, Ly) in [(3, 3), (3, 6), (6, 6)]
            λ = QAtlas.fetch(
                QAtlas.Triangular(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=1.0
            )
            @test count(x -> abs(x - 3.0) < 1e-10, λ) >= 2
            @test maximum(λ) ≈ 3.0 atol = 1e-10
        end
    end

    @testset "Frustration: spectrum NOT symmetric about zero" begin
        # Unlike honeycomb / Lieb (bipartite), the triangular band is
        # asymmetric: min = -6, max = +3 (not ±6).
        λ = QAtlas.fetch(QAtlas.Triangular(), TightBindingSpectrum(); Lx=6, Ly=6, t=1.0)
        @test minimum(λ) ≈ -6.0 atol = 1e-10
        @test maximum(λ) ≈ +3.0 atol = 1e-10
        @test !(sort(λ) ≈ sort(-λ))  # not chiral-symmetric
    end

    @testset "t scaling: λ(t) = t · λ(1)" begin
        λ1 = QAtlas.fetch(QAtlas.Triangular(), TightBindingSpectrum(); Lx=4, Ly=4, t=1.0)
        for t in (0.5, 2.0, 3.7)
            λt = QAtlas.fetch(QAtlas.Triangular(), TightBindingSpectrum(); Lx=4, Ly=4, t=t)
            @test λt ≈ t .* λ1 rtol = 1e-12
        end
    end
end
