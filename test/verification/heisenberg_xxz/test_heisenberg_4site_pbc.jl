# ─────────────────────────────────────────────────────────────────────────────
# Verification: 4-site spin-1/2 Heisenberg ring (PBC)
#
# Build H = J Σ_{i mod 4} S_i · S_{i+1} via Lattice2D's 4-site PBC
# square chain (build_lattice(Square, 4, 1; PeriodicAxis in x, OpenAxis
# in y)) and compare the full 2^4 = 16-state spectrum against the
# hardcoded exact result in src/models/quantum/Heisenberg/Heisenberg.jl.
#
# Two independent construction paths cross-check each other:
#   src:  hardcoded {-2J, -J×3, 0×7, +J×5} (from textbook / Bethe ansatz)
#   test: Lattice2D bonds(lat) → build_spinhalf_heisenberg → eigvals
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

@testset "Heisenberg 4-site PBC — ED vs exact spectrum" begin
    # 4-site PBC ring via Lattice2D: PBC in x, OBC in y (drops the
    # (0,1) connection since Ly = 1).
    lat = build_lattice(
        Square, 4, 1; boundary=LatticeBoundary((PeriodicAxis(), OpenAxis()))
    )
    @test num_sites(lat) == 4
    @test length(collect(bonds(lat))) == 4

    @testset "Full spectrum for various J" begin
        for J in (1.0, 0.5, 2.0, -1.0)
            H = build_spinhalf_heisenberg(lat, J)
            λ_ed = sort(eigvals(Symmetric(H)))
            λ_exact = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=4, J=J, bc=:PBC)
            @test length(λ_ed) == 16
            @test length(λ_exact) == 16
            @test λ_ed ≈ λ_exact atol = 1e-12
        end
    end

    @testset "Ground state E_0 = -2J (unique singlet)" begin
        for J in (1.0, 0.5, 3.7)
            λ = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=4, J=J, bc=:PBC)
            @test λ[1] ≈ -2J atol = 1e-12
            @test λ[2] > λ[1] + 1e-10
        end
    end

    @testset "Degeneracy structure" begin
        λ = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=4, J=1.0, bc=:PBC)
        @test count(x -> abs(x + 2.0) < 1e-10, λ) == 1   # singlet
        @test count(x -> abs(x + 1.0) < 1e-10, λ) == 3   # triplet
        @test count(x -> abs(x) < 1e-10, λ) == 7          # singlet + 2 triplets
        @test count(x -> abs(x - 1.0) < 1e-10, λ) == 5   # quintet
    end

    @testset "Structural: tr H = 0 and Hermitian" begin
        H = build_spinhalf_heisenberg(lat, 1.0)
        @test tr(H) ≈ 0 atol = 1e-14
        @test H ≈ H' atol = 1e-14
    end

    @testset "Quintet at E = J (ferromagnetic sector)" begin
        # All spins up: ⟨↑↑↑↑|H|↑↑↑↑⟩ = J·4·(1/4) = J
        # (4 bonds, each Sz·Sz = 1/4, ladder terms vanish)
        H = build_spinhalf_heisenberg(lat, 1.0)
        all_up = zeros(16)
        all_up[1] = 1.0   # |↑↑↑↑⟩ in kron convention
        @test (H * all_up) ≈ (1.0 * all_up) atol = 1e-14
    end

    @testset "E_0 per site approaches Bethe-ansatz limit" begin
        # For 4-site PBC: E₀ / N = -2J / 4 = -0.5J
        # Bethe ansatz (N → ∞): e₀ = (1/4 - ln 2)J ≈ -0.4431J
        # 4-site overshoot is expected (finite-size effect).
        λ = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=4, J=1.0, bc=:PBC)
        e0_per_site = λ[1] / 4
        @test e0_per_site < (0.25 - log(2))  # below Bethe ansatz
        @test e0_per_site > -0.6              # not too far off
    end
end
