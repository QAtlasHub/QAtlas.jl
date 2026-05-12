# ─────────────────────────────────────────────────────────────────────────────
# Verification: 2-site spin-1/2 Heisenberg dimer — full spectrum
#
# Compare Kronecker-embed ED of the real-space Hamiltonian
#     H = J S_1 · S_2
# built from `Lattice2D`'s `bonds(lat)` to the analytical spectrum
#     { -3J/4, J/4, J/4, J/4 }  (singlet + three-fold triplet)
# recorded in `src/models/quantum/Heisenberg/Heisenberg.jl`.
#
# Also verifies the singlet–triplet gap Δ = J and that the low-lying
# eigenvalue signs flip correctly under J → −J (ferromagnetic flip).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

@testset "Heisenberg dimer (N=2) — ED vs exact spectrum" begin
    # A 2-site OBC chain: build_lattice(Square, 2, 1) with OBC in y drops
    # the (0, 1) connection; only the single (1, 0) bond survives.
    lat = build_lattice(Square, 2, 1; boundary=OpenAxis())
    @test num_sites(lat) == 2
    @test length(collect(bonds(lat))) == 1

    @testset "Full spectrum for various J" begin
        for J in (1.0, 0.5, 2.0, -1.0)
            H = build_spinhalf_heisenberg(lat, J)
            λ_ed = sort(eigvals(Symmetric(H)))
            λ_exact = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=2, J=J)
            @test λ_ed ≈ λ_exact atol = 1e-12
        end
    end

    @testset "Singlet-triplet gap Δ = J" begin
        for J in (1.0, 0.5, 3.7)
            λ = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=2, J=J)
            # Sorted: J>0 → [-3J/4, J/4, J/4, J/4]; gap = λ[2] - λ[1]
            gap = λ[2] - λ[1]
            @test gap ≈ J rtol = 1e-12
        end
    end

    @testset "Structural: tr H = 0 and Hermitian" begin
        H = build_spinhalf_heisenberg(lat, 1.0)
        @test tr(H) ≈ 0 atol = 1e-14
        @test H ≈ H' atol = 1e-14
    end

    @testset "Singlet is the ground state for J > 0" begin
        # Eigenvector at −3J/4 should be (|↑↓⟩ − |↓↑⟩)/√2 (up to phase).
        # Basis ordering (kron with site 1 outermost):
        #   1 = |↑↑⟩, 2 = |↑↓⟩, 3 = |↓↑⟩, 4 = |↓↓⟩
        H = build_spinhalf_heisenberg(lat, 1.0)
        F = eigen(Symmetric(H))
        @test F.values[1] ≈ -0.75 atol = 1e-12
        ψ0 = F.vectors[:, 1]
        # |c(|↑↑⟩)|² = |c(|↓↓⟩)|² = 0
        @test abs2(ψ0[1]) < 1e-20
        @test abs2(ψ0[4]) < 1e-20
        # ψ0 antisymmetric between |↑↓⟩ and |↓↑⟩: c[2] = -c[3] (up to sign).
        @test ψ0[2] ≈ -ψ0[3] atol = 1e-12
        # Singlet normalization: |c[2]| = |c[3]| = 1/√2.
        @test abs(ψ0[2]) ≈ 1 / sqrt(2) atol = 1e-12
    end

    @testset "Ferromagnetic flip J → -J inverts ordering" begin
        # For J < 0, the triplet becomes the ground state (E_t = J/4 < 0).
        λ_afm = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=2, J=1.0)
        λ_fm = QAtlas.fetch(Heisenberg1D(), ExactSpectrum(); N=2, J=-1.0)
        @test λ_afm ≈ -reverse(λ_fm) atol = 1e-12
    end
end
