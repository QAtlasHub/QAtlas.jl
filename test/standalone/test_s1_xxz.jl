using QAtlas, Test, LinearAlgebra

# Spin-1 XXZ chain — small-N dense-ED reference for ThermalMPS validation.
# Hilbert space 3^N (capped at N ≤ 8).
#
# Convention: spin-1 operators (eigenvalues ±1, 0).
#   H = J Σᵢ [ Sˣ Sˣ + Sʸ Sʸ + Δ Sᶻ Sᶻ ]
# At Δ = 1 this coincides with S1Heisenberg1D; at Δ = 0 it is the spin-1 XY model.

@testset "S1XXZ1D — Δ = 1 reproduces S1Heisenberg1D exactly" begin
    # At Δ = 1, S1XXZ1D Hamiltonian equals S1Heisenberg1D's bond-by-bond.
    for (J, N) in ((1.0, 3), (1.3, 4), (0.7, 4))
        H_xxz = QAtlas._s1_xxz_hamiltonian_matrix(S1XXZ1D(; J=J, Δ=1.0), N)
        H_hei = QAtlas._s1_heisenberg_hamiltonian_matrix(S1Heisenberg1D(; J=J), N)
        @test maximum(abs.(H_xxz - H_hei)) < 1e-12
    end
end

@testset "S1XXZ1D — Δ = 1 thermal energy matches S1Heisenberg1D fetch" begin
    for (J, N, β) in ((1.0, 3, 1.0), (1.3, 4, 0.7))
        E_xxz = QAtlas.fetch(S1XXZ1D(; J=J, Δ=1.0), Energy(), OBC(N); beta=β)
        E_hei = QAtlas.fetch(S1Heisenberg1D(; J=J), Energy(), OBC(N); beta=β)
        @test E_xxz ≈ E_hei rtol = 1e-12
    end
end

@testset "S1XXZ1D — Tr(H) = 0 (β → 0 ⇒ ⟨H⟩ = 0)" begin
    # Each Sᵅ is traceless and Sᵅ ⊗ Sᵅ has Tr·Tr = 0.
    for Δ in (-0.5, 0.0, 0.5, 1.0), N in (2, 3, 4)
        E0 = QAtlas.fetch(S1XXZ1D(; J=1.0, Δ=Δ), Energy(), OBC(N); beta=0.0)
        @test abs(E0) < 1e-12
    end
end

@testset "S1XXZ1D — β → ∞ recovers OBC ground state" begin
    for (J, Δ, N) in ((1.0, 0.5, 4), (1.0, -0.3, 4), (1.3, 1.0, 3))
        m = S1XXZ1D(; J=J, Δ=Δ)
        H = QAtlas._s1_xxz_hamiltonian_matrix(m, N)
        evals = sort(real.(eigvals(Hermitian(H))))
        E_gs = evals[1]
        gap = evals[2] - evals[1]
        β = 50.0
        E_th = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        tol = max(1e-10, 10 * gap * exp(-β * gap))
        @test abs(E_th - E_gs) ≤ tol
    end
end

@testset "S1XXZ1D — direct ED cross-check at small N (β finite)" begin
    for (J, Δ, N, β) in ((1.0, 0.5, 3, 1.0), (1.3, -0.5, 4, 0.7), (0.6, 0.0, 4, 2.5))
        m = S1XXZ1D(; J=J, Δ=Δ)
        H = QAtlas._s1_xxz_hamiltonian_matrix(m, N)
        evals = real.(eigvals(Hermitian(H)))
        emin = minimum(evals)
        ws = exp.(-β .* (evals .- emin))
        E_direct = sum(evals .* ws) / sum(ws)
        E_fetch = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        @test E_fetch ≈ E_direct rtol = 1e-12
    end
end

@testset "S1XXZ1D — ε = f + T·s identity (per-site)" begin
    for (J, Δ, N, β) in ((1.0, 0.5, 3, 0.5), (1.0, -0.3, 4, 1.0), (1.3, 1.0, 4, 1.5))
        m = S1XXZ1D(; J=J, Δ=Δ)
        E_total = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        f = QAtlas.fetch(m, FreeEnergy(), OBC(N); beta=β)
        s = QAtlas.fetch(m, ThermalEntropy(), OBC(N); beta=β)
        ε = E_total / N
        @test ε ≈ f + s / β rtol = 1e-12
    end
end

@testset "S1XXZ1D — SpecificHeat = β²·Var(H)/N" begin
    for (J, Δ, N, β) in ((1.0, 0.5, 3, 1.0), (1.3, -0.5, 4, 0.5))
        m = S1XXZ1D(; J=J, Δ=Δ)
        H = QAtlas._s1_xxz_hamiltonian_matrix(m, N)
        evals = real.(eigvals(Hermitian(H)))
        emin = minimum(evals)
        ws = exp.(-β .* (evals .- emin))
        Z = sum(ws)
        E1 = sum(evals .* ws) / Z
        E2 = sum((evals .^ 2) .* ws) / Z
        c_direct = β^2 * (E2 - E1^2) / N
        c_fetch = QAtlas.fetch(m, SpecificHeat(), OBC(N); beta=β)
        @test c_fetch ≈ c_direct rtol = 1e-12
    end
end

@testset "S1XXZ1D — ⟨Sᵅ⟩ = 0 by SU(2) / U(1) symmetry" begin
    # H is U(1)-symmetric in Sᶻ for any Δ (preserves Σ Sᶻ); at Δ = 1 it is
    # SU(2)-symmetric.  Without symmetry-breaking field the bulk
    # magnetisations vanish at any β.
    for (Δ, N, β) in ((0.5, 3, 1.0), (-0.3, 4, 0.7), (1.0, 4, 0.5))
        m = S1XXZ1D(; J=1.0, Δ=Δ)
        @test abs(QAtlas.fetch(m, MagnetizationX(), OBC(N); beta=β)) < 1e-12
        @test abs(QAtlas.fetch(m, MagnetizationY(), OBC(N); beta=β)) < 1e-12
        @test abs(QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β)) < 1e-12
    end
end

@testset "S1XXZ1D — Δ = 1 correlator equals S1Heisenberg1D correlator" begin
    # Sanity: <S^α S^α> at Δ=1 must agree with the dedicated Heisenberg path.
    m_xxz = S1XXZ1D(; J=1.0, Δ=1.0)
    m_hei = S1Heisenberg1D(; J=1.0)
    N = 3
    β = 1.0
    for axis in (XXCorrelation, YYCorrelation, ZZCorrelation), i in 1:N, j in 1:N
        c1 = QAtlas.fetch(m_xxz, axis{:static}(), OBC(N); beta=β, i=i, j=j)
        c2 = QAtlas.fetch(m_hei, axis{:static}(), OBC(N); beta=β, i=i, j=j)
        @test c1 ≈ c2 atol = 1e-12
    end
end

@testset "S1XXZ1D — connected = static when one operator has zero mean" begin
    # SU(2)/U(1) gives ⟨Sᵅ⟩ = 0, so connected XX/YY/ZZ correlators equal
    # the static ones in the symmetric ensemble.
    m = S1XXZ1D(; J=1.0, Δ=0.5)
    N = 3
    β = 1.0
    for axis in (XXCorrelation, YYCorrelation, ZZCorrelation), i in 1:N, j in 1:N
        c_st = QAtlas.fetch(m, axis{:static}(), OBC(N); beta=β, i=i, j=j)
        c_co = QAtlas.fetch(m, axis{:connected}(), OBC(N); beta=β, i=i, j=j)
        @test c_st ≈ c_co atol = 1e-12
    end
end

@testset "S1XXZ1D — MassGap matches direct E₁ - E₀" begin
    for (J, Δ, N) in ((1.0, 0.5, 3), (1.0, -0.3, 4), (1.3, 1.0, 4))
        m = S1XXZ1D(; J=J, Δ=Δ)
        H = QAtlas._s1_xxz_hamiltonian_matrix(m, N)
        evals = sort(real.(eigvals(Hermitian(H))))
        Δ_direct = evals[2] - evals[1]
        Δ_fetch = QAtlas.fetch(m, MassGap(), OBC(N))
        @test Δ_fetch ≈ Δ_direct rtol = 1e-12
    end
end

@testset "S1XXZ1D — ExactSpectrum returns sorted 3^N eigenvalues" begin
    for (Δ, N) in ((0.5, 3), (-0.3, 4))
        m = S1XXZ1D(; J=1.0, Δ=Δ)
        spec = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
        @test length(spec) == 3^N
        @test issorted(spec)
    end
end

@testset "S1XXZ1D — N bounds enforced" begin
    @test_throws ArgumentError QAtlas._s1_xxz_hamiltonian_matrix(S1XXZ1D(; J=1.0, Δ=0.5), 1)
    @test_throws ArgumentError QAtlas._s1_xxz_hamiltonian_matrix(S1XXZ1D(; J=1.0, Δ=0.5), 9)
end
