using QAtlas, Test, LinearAlgebra

# Cross-check the dense-ED thermal Energy on XXZ1D OBC against an
# independent direct ED at small N and against sanity limits.

@testset "XXZ1D thermal OBC: β → 0 gives zero energy (Tr(H) = 0 on spin-1/2)" begin
    # Every bond term σᵅ⊗σᵅ is traceless, and σᶻ at any site has Tr = 0,
    # so the XXZ OBC Hamiltonian satisfies Tr(H) = 0 exactly.  At β = 0
    # this means `⟨H⟩ = Tr(H) / 2^N = 0`.
    for Δ in (-0.5, 0.0, 0.7, 1.0), N in (4, 5, 6)
        E_total = QAtlas.fetch(XXZ1D(; J=1.0, Δ=Δ), Energy(), OBC(N); beta=0.0)
        @test abs(E_total) < 1e-12
    end
end

@testset "XXZ1D thermal OBC: β → ∞ reproduces the OBC GS energy (by diagonalisation)" begin
    # Compute the true OBC GS energy by diagonalising the same matrix
    # and verify that a large-β thermal average collapses onto it up to
    # `exp(-β · Δ_gap)`.  Δ_gap is determined empirically from the
    # spectrum.
    for (J, Δ, N) in ((1.0, 1.0, 5), (1.3, 0.5, 5))
        model = XXZ1D(; J=J, Δ=Δ)
        H = QAtlas._xxz1d_hamiltonian_matrix(model, N)
        evals = sort(real.(eigvals(Hermitian(H))))
        E_gs = evals[1]
        gap = evals[2] - evals[1]
        β = 50.0
        E_thermal = QAtlas.fetch(model, Energy(), OBC(N); beta=β)  # total energy
        # At β large, correction ≈ gap · exp(-β gap), bounded very loosely.
        tol = max(1e-10, 10 * gap * exp(-β * gap))
        @test abs(E_thermal - E_gs) ≤ tol
    end
end

@testset "XXZ1D thermal OBC: matches direct ED with Pauli bond matrix (N=2, N=3)" begin
    # Direct independent construction for very small N so the kron chain
    # is unambiguous.
    σx = ComplexF64[0 1; 1 0]
    σy = ComplexF64[0 -im; im 0]
    σz = ComplexF64[1 0; 0 -1]
    I2 = ComplexF64[1 0; 0 1]

    for (J, Δ, β) in ((1.3, 0.7, 2.0), (0.8, -0.5, 0.7))
        # N = 2: only one bond, H is just the bond block.
        bond_22 = (J / 4) * (kron(σx, σx) + kron(σy, σy) + Δ * kron(σz, σz))
        H2 = Hermitian(bond_22)
        evals2 = real.(eigvals(H2))
        emin2 = minimum(evals2)
        ws2 = exp.(-β .* (evals2 .- emin2))
        E2_direct = real(sum(evals2 .* ws2) / sum(ws2))
        E2_fetch = QAtlas.fetch(XXZ1D(; J=J, Δ=Δ), Energy(), OBC(2); beta=β)
        @test E2_fetch ≈ E2_direct rtol = 1e-10

        # N = 3: two bonds, H = bond(1,2)⊗I + I⊗bond(2,3).
        bond12 = kron(bond_22, I2)
        bond23 = kron(I2, bond_22)
        H3 = Hermitian(bond12 + bond23)
        evals3 = real.(eigvals(H3))
        emin3 = minimum(evals3)
        ws3 = exp.(-β .* (evals3 .- emin3))
        E3_direct = real(sum(evals3 .* ws3) / sum(ws3))
        E3_fetch = QAtlas.fetch(XXZ1D(; J=J, Δ=Δ), Energy(), OBC(3); beta=β)
        @test E3_fetch ≈ E3_direct rtol = 1e-10
    end
end

@testset "XXZ1D thermal OBC: monotone cooling (E(β₁) ≥ E(β₂) for β₁ < β₂)" begin
    N = 5
    model = XXZ1D(; J=1.0, Δ=0.5)
    E0 = QAtlas.fetch(model, Energy(), OBC(N); beta=0.0)
    E1 = QAtlas.fetch(model, Energy(), OBC(N); beta=0.5)
    E2 = QAtlas.fetch(model, Energy(), OBC(N); beta=1.0)
    E5 = QAtlas.fetch(model, Energy(), OBC(N); beta=5.0)
    @test E0 ≥ E1 ≥ E2 ≥ E5
    @test isapprox(E0, 0.0; atol=1e-12)
end

@testset "Heisenberg1D thermal OBC: matches XXZ1D(Δ=1) thermal OBC" begin
    # The delegator should produce the same number as XXZ1D(Δ=1) at
    # matching J.
    for β in (0.3, 1.0, 2.5), N in (4, 6)
        E_heis = QAtlas.fetch(Heisenberg1D(), Energy(), OBC(N); beta=β, J=1.5)
        E_xxz = QAtlas.fetch(XXZ1D(; J=1.5, Δ=1.0), Energy(), OBC(N); beta=β)
        @test E_heis ≈ E_xxz rtol = 1e-12
    end
end

@testset "dense-ED cap guards against runaway N" begin
    @test_throws ArgumentError QAtlas._pauli_string(
        QAtlas._MAX_ED_SITES + 1, 1 => QAtlas._σx
    )
end
