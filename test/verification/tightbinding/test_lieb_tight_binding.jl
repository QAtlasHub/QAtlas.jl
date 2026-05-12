# ─────────────────────────────────────────────────────────────────────────────
# Verification: nearest-neighbor tight-binding on the Lieb lattice
#
# Cross-validate QAtlas's closed-form Bloch spectrum (3-band per k)
# against direct diagonalization of the real-space tight-binding
# Hamiltonian built from `Lattice2D.bonds(lat)`.
#
# Test sizes: 2×2 through 4×4 Lieb PBC.
#
# Additional structural checks:
#   * tr H = 0                              (NN TB has no on-site energy)
#   * Spectrum symmetric about zero         (bipartite A vs {B, C})
#   * Lx·Ly zero modes (generic flat band)  — plus 2 extra at the
#     M-point if both Lx and Ly are even    (three-fold band touching)
#   * Ground state at -2√2·t reached at Γ   (cos²(0)+cos²(0) = 2)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

const T_HOP = 1.0

expected_zero_count(Lx, Ly) = Lx * Ly + (iseven(Lx) && iseven(Ly) ? 2 : 0)

@testset "Lieb TB — real-space vs Bloch closed form" begin
    for (Lx, Ly) in [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]
        lat = build_lattice(Lieb, Lx, Ly)

        @testset "$(Lx)×$(Ly) Lieb PBC" begin
            H = build_tight_binding(lat, T_HOP)
            λ_real = sort(eigvals(Symmetric(H)))
            λ_bloch = QAtlas.fetch(
                QAtlas.Lieb(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=T_HOP
            )

            @test length(λ_real) == 3 * Lx * Ly
            @test length(λ_bloch) == 3 * Lx * Ly
            @test λ_real ≈ λ_bloch atol = 1e-10

            # Structural: tr H = 0
            @test tr(H) ≈ 0 atol = 1e-10

            # Bipartite spectrum: {λ} is symmetric about zero
            @test sort(λ_real) ≈ sort(-λ_real) atol = 1e-10

            # Flat band at zero: Lx·Ly states, plus 2 extra at M-point if
            # both Lx and Ly are even.
            zeros = count(x -> abs(x) < 1e-10, λ_real)
            @test zeros == expected_zero_count(Lx, Ly)

            # Ground state: -2√2·t at the Γ-point
            @test λ_real[1] ≈ -2 * sqrt(2) atol = 1e-10
        end
    end

    @testset "2×2 hand-computable spectrum" begin
        # Allowed (θ₁, θ₂) ∈ {0, π}²:
        #   (0, 0) : E = ±2·√(1 + 1)·t = ±2√2·t, plus 0   → {-2√2, 0, 2√2}
        #   (π, 0) : E = ±2·√(0 + 1)·t = ±2t, plus 0      → {-2, 0, 2}
        #   (0, π) : identical to (π, 0) by symmetry      → {-2, 0, 2}
        #   (π, π) : E = 0, plus 0 → all three zero       → {0, 0, 0}
        # Sorted full spectrum: [-2√2, -2, -2, 0×6, 2, 2, 2√2]
        λ = QAtlas.fetch(QAtlas.Lieb(), TightBindingSpectrum(); Lx=2, Ly=2, t=1.0)
        expected = sort([
            -2 * sqrt(2), -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2 * sqrt(2)
        ])
        @test λ ≈ expected atol = 1e-12
    end

    @testset "3×3 Lieb — no M-point zero modes" begin
        # M-point (θ=π each) requires both Lx, Ly even. For 3×3 the
        # spectrum has exactly Lx·Ly = 9 zero eigenvalues (flat band
        # only, no band-touching bonus).
        λ = QAtlas.fetch(QAtlas.Lieb(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
        @test count(x -> abs(x) < 1e-10, λ) == 9
    end

    @testset "t scaling: λ(t) = t · λ(1)" begin
        λ1 = QAtlas.fetch(QAtlas.Lieb(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)
        for t in (0.5, 2.0, 3.7)
            λt = QAtlas.fetch(QAtlas.Lieb(), TightBindingSpectrum(); Lx=3, Ly=3, t=t)
            @test λt ≈ t .* λ1 rtol = 1e-12
        end
    end
end
