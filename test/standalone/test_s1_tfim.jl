using QAtlas, Test, LinearAlgebra

# Spin-1 transverse-field Ising chain — small-N dense-ED reference for
# ThermalMPS validation.  Hilbert space 3^N (capped at N ≤ 8).
#
# Conventions (from `_S1_x`, `_S1_z`):
#   spin-1 operators with eigenvalues ±1, 0
#   H = -J Σᵢ Sᶻᵢ Sᶻᵢ₊₁ - h Σᵢ Sˣᵢ

@testset "S1TFIM — N=2 spectrum at h = 0 (classical Ising bond)" begin
    # At h = 0, H = -J Sᶻ⊗Sᶻ.  Sᶻ has eigenvalues {+1, 0, -1}, so Sᶻ⊗Sᶻ
    # has eigenvalues {+1 (×2: 1·1, -1·-1), 0 (×5), -1 (×2: 1·-1, -1·1)}.
    # H = -J · (those) gives [-J ×2, 0 ×5, +J ×2].
    for J in (0.7, 1.0, 1.3)
        H = QAtlas._s1_tfim_hamiltonian_matrix(S1TFIM(; J=J, h=0.0), 2)
        evals = sort(real.(eigvals(Hermitian(H))))
        expected = sort([fill(-J, 2); fill(0.0, 5); fill(J, 2)])
        @test evals ≈ expected atol = 1e-12
    end
end

@testset "S1TFIM — N=1 transverse-field-only single site at J → 0" begin
    # Single-bond chain (N=2) at J = 0 is just -h (Sˣ ⊗ I + I ⊗ Sˣ).
    # Sˣ has eigenvalues {-1, 0, +1}, so the spectrum is the pairwise sum:
    # {-2, -1, -1, 0, 0, 0, +1, +1, +2} times -h.
    for h in (0.5, 1.0, 1.7)
        H = QAtlas._s1_tfim_hamiltonian_matrix(S1TFIM(; J=0.0, h=h), 2)
        evals = sort(real.(eigvals(Hermitian(H))))
        expected = sort(-h .* [-2.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0])
        @test evals ≈ expected atol = 1e-10
    end
end

@testset "S1TFIM — β → ∞ recovers OBC ground state" begin
    for (J, h, N) in ((1.0, 0.5, 3), (1.0, 1.0, 4), (0.7, 0.3, 4))
        m = S1TFIM(; J=J, h=h)
        H = QAtlas._s1_tfim_hamiltonian_matrix(m, N)
        evals = sort(real.(eigvals(Hermitian(H))))
        E_gs = evals[1]
        gap = evals[2] - evals[1]
        β = 50.0
        E_th = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        tol = max(1e-10, 10 * gap * exp(-β * gap))
        @test abs(E_th - E_gs) ≤ tol
    end
end

@testset "S1TFIM — direct ED cross-check at small N (β finite)" begin
    for (J, h, N, β) in ((1.0, 0.5, 3, 1.0), (1.3, 0.7, 4, 0.7), (0.6, 1.5, 4, 2.5))
        m = S1TFIM(; J=J, h=h)
        H = QAtlas._s1_tfim_hamiltonian_matrix(m, N)
        evals = real.(eigvals(Hermitian(H)))
        emin = minimum(evals)
        ws = exp.(-β .* (evals .- emin))
        E_direct = sum(evals .* ws) / sum(ws)
        E_fetch = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        @test E_fetch ≈ E_direct rtol = 1e-12
    end
end

@testset "S1TFIM — ε = f + T·s identity (per-site)" begin
    for (J, h, N, β) in ((1.0, 0.5, 3, 0.5), (1.0, 1.0, 4, 1.0), (0.7, 0.3, 4, 1.5))
        m = S1TFIM(; J=J, h=h)
        E_total = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        f = QAtlas.fetch(m, FreeEnergy(), OBC(N); beta=β)
        s = QAtlas.fetch(m, ThermalEntropy(), OBC(N); beta=β)
        ε = E_total / N            # per-site
        @test ε ≈ f + s / β rtol = 1e-12
    end
end

@testset "S1TFIM — SpecificHeat = β² · Var(H) / N" begin
    for (J, h, N, β) in ((1.0, 0.5, 3, 1.0), (1.0, 1.0, 4, 0.5))
        m = S1TFIM(; J=J, h=h)
        H = QAtlas._s1_tfim_hamiltonian_matrix(m, N)
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

@testset "S1TFIM — ⟨Sᶻ⟩ = 0 by Z₂ symmetry (h ≠ 0)" begin
    # Z₂ flip P: Sˣ → Sˣ, Sᶻ → -Sᶻ commutes with H, so ⟨Sᶻ⟩_β = 0 in
    # the symmetric (finite-β, finite-N) phase.
    for (J, h, N, β) in ((1.0, 0.5, 3, 1.0), (1.0, 1.0, 4, 0.5), (0.7, 1.5, 3, 2.0))
        m = S1TFIM(; J=J, h=h)
        mz = QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β)
        @test abs(mz) < 1e-12
    end
end

@testset "S1TFIM — ⟨(Sᶻ)²⟩ = ⟨Sᶻ_i Sᶻ_i⟩ self-correlator" begin
    # ZZCorrelation{:static}(i, i) should equal Tr(ρ · _spin1_string(N, i => Sᶻ²))
    # which is bounded between 0 and 1 (eigenvalues of (Sᶻ)² are 0 or 1).
    m = S1TFIM(; J=1.0, h=0.5)
    N = 3
    β = 1.0
    for i in 1:N
        c = QAtlas.fetch(m, ZZCorrelation{:static}(), OBC(N); beta=β, i=i, j=i)
        @test 0.0 ≤ c ≤ 1.0 + 1e-12
    end
end

@testset "S1TFIM — ZZCorrelation connected = 0 at h = 0 (classical limit)" begin
    # At h = 0 the eigenstates are product states in the Sᶻ basis.  At
    # any β the thermal state is diagonal in Sᶻ so ⟨Sᶻ⟩ = 0 by Z₂ and
    # ⟨Sᶻ_i Sᶻ_j⟩ - ⟨Sᶻ_i⟩⟨Sᶻ_j⟩ = ⟨Sᶻ_i Sᶻ_j⟩.
    m = S1TFIM(; J=1.0, h=0.0)
    N = 3
    β = 1.0
    for i in 1:N, j in 1:N
        c_static = QAtlas.fetch(m, ZZCorrelation{:static}(), OBC(N); beta=β, i=i, j=j)
        c_conn   = QAtlas.fetch(m, ZZCorrelation{:connected}(), OBC(N); beta=β, i=i, j=j)
        @test c_conn ≈ c_static atol = 1e-12  # ⟨Sᶻ⟩ = 0 ⇒ connected = static
    end
end

@testset "S1TFIM — MassGap matches direct E₁ - E₀" begin
    for (J, h, N) in ((1.0, 0.5, 3), (1.0, 1.0, 4))
        m = S1TFIM(; J=J, h=h)
        H = QAtlas._s1_tfim_hamiltonian_matrix(m, N)
        evals = sort(real.(eigvals(Hermitian(H))))
        Δ_direct = evals[2] - evals[1]
        Δ_fetch = QAtlas.fetch(m, MassGap(), OBC(N))
        @test Δ_fetch ≈ Δ_direct rtol = 1e-12
    end
end

@testset "S1TFIM — ExactSpectrum returns sorted 3^N eigenvalues" begin
    for (J, h, N) in ((1.0, 0.5, 2), (1.0, 1.0, 3), (0.7, 0.3, 4))
        m = S1TFIM(; J=J, h=h)
        spec = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
        @test length(spec) == 3^N
        @test issorted(spec)
        # First eigenvalue equals the GS energy.
        H = QAtlas._s1_tfim_hamiltonian_matrix(m, N)
        @test spec[1] ≈ minimum(real.(eigvals(Hermitian(H)))) atol = 1e-12
    end
end

@testset "S1TFIM — N bounds enforced" begin
    @test_throws ArgumentError QAtlas._s1_tfim_hamiltonian_matrix(S1TFIM(; J=1.0, h=0.5), 1)
    @test_throws ArgumentError QAtlas._s1_tfim_hamiltonian_matrix(S1TFIM(; J=1.0, h=0.5), 9)
end