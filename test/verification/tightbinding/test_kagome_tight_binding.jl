# ─────────────────────────────────────────────────────────────────────────────
# Verification: nearest-neighbor tight-binding on the kagome lattice
#
# Cross-validate QAtlas's closed-form Bloch spectrum (3×3 per k) against
# direct diagonalization of the real-space tight-binding Hamiltonian
# built from `Lattice2D.bonds(lat)`.
#
# Test sizes: 2×2 through 4×4 kagome PBC.
#
# Additional structural checks:
#   * tr H = 0 (NN TB has no on-site energy)
#   * Exactly Lx·Ly (+1 from Γ-touch) eigenvalues at +2t  → flat-band
#     existence
#   * Bonding ground state E_gs = −4t, reached at the Γ-point
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

const T_HOP = 1.0

@testset "Kagome TB — real-space vs Bloch closed form" begin
    for (Lx, Ly) in [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]
        lat = build_lattice(Kagome, Lx, Ly)

        @testset "$(Lx)×$(Ly) kagome PBC" begin
            H = build_tight_binding(lat, T_HOP)
            λ_real = sort(eigvals(Symmetric(H)))
            λ_bloch = QAtlas.fetch(
                QAtlas.Kagome(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=T_HOP
            )

            @test length(λ_real) == 3 * Lx * Ly
            @test length(λ_bloch) == 3 * Lx * Ly
            @test λ_real ≈ λ_bloch atol = 1e-10

            # Structural: trH = 0
            @test tr(H) ≈ 0 atol = 1e-10

            # Flat band at +2t: Lx·Ly states from the flat band, plus
            # one extra from the Γ-point band touching.
            flat = count(x -> abs(x - 2.0) < 1e-10, λ_real)
            @test flat == Lx * Ly + 1

            # Ground state: -4t at the Γ-point (unique)
            @test λ_real[1] ≈ -4.0 atol = 1e-10
            @test count(x -> abs(x - (-4.0)) < 1e-10, λ_real) == 1
        end
    end

    @testset "2×2 hand-computable spectrum" begin
        # Allowed k-points: (θ₁, θ₂) ∈ {0, π} × {0, π}
        #   (0, 0)   : eigs of -2·[[0,1,1],[1,0,1],[1,1,0]]        = {-4, +2, +2}
        #   (π, 0)   : eigs of -2·[[0,0,1],[0,0, cos(-π/2)=0],...] = {-2, 0, +2}
        #   (0, π)   :                                              = {-2, 0, +2}
        #   (π, π)   :                                              = {-2, 0, +2}
        # Sorted full spectrum: [-4, -2, -2, -2, 0, 0, 0, 2, 2, 2, 2, 2]
        λ = QAtlas.fetch(QAtlas.Kagome(), TightBindingSpectrum(); Lx=2, Ly=2, t=1.0)
        @test λ ≈ [-4.0, -2.0, -2.0, -2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0] atol =
            1e-12
    end

    @testset "t scaling: λ(t) = t · λ(1)" begin
        λ1 = QAtlas.fetch(QAtlas.Kagome(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
        for t in (0.5, 2.0, 3.7)
            λt = QAtlas.fetch(QAtlas.Kagome(), TightBindingSpectrum(); Lx=3, Ly=3, t=t)
            @test λt ≈ t .* λ1 rtol = 1e-12
        end
    end
end
