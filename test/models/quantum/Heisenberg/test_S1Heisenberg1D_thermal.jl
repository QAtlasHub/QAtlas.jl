using QAtlas, Test, LinearAlgebra

# Spin-1 Heisenberg (Haldane chain) — small-N dense-ED reference for
# ThermalMPS validation.  Bond eigenvalues from S_tot ∈ {0,1,2}:
#   singlet = -2J,  triplet = -J×3,  quintet = +J×5  (per bond).
#
# Total Hilbert space 3^N (capped at N ≤ 8 by `_MAX_ED_SITES_S1`).

@testset "S1Heisenberg1D — N=2 dimer spectrum (analytic)" begin
    # Two spin-1 → 9 states: [-2J, -J×3, +J×5].
    for J in (0.7, 1.0, 1.3)
        H = QAtlas._s1_heisenberg_hamiltonian_matrix(S1Heisenberg1D(; J=J), 2)
        evals = sort(real.(eigvals(Hermitian(H))))
        expected = sort([-2J; fill(-J, 3); fill(J, 5)])
        @test evals ≈ expected atol = 1e-12
    end
end

@testset "S1Heisenberg1D — Tr(H) = 0  (β → 0 gives ⟨H⟩ = 0)" begin
    # Sᵅ ⊗ Sᵅ has trace Tr(Sᵅ)·Tr(Sᵅ) = 0 for every α (each Sᵅ is traceless).
    for J in (0.5, 1.0), N in (2, 3, 4)
        E0 = QAtlas.fetch(S1Heisenberg1D(; J=J), Energy(), OBC(N); beta=0.0)
        @test abs(E0) < 1e-12
    end
end

@testset "S1Heisenberg1D — β → ∞ recovers OBC ground state" begin
    # Compute true ground from the same Hamiltonian and check the large-β
    # thermal mean collapses onto E_gs up to exp(-β · gap).
    for (J, N) in ((1.0, 3), (1.0, 4), (0.7, 4))
        m = S1Heisenberg1D(; J=J)
        H = QAtlas._s1_heisenberg_hamiltonian_matrix(m, N)
        evals = sort(real.(eigvals(Hermitian(H))))
        E_gs = evals[1]
        gap = evals[2] - evals[1]
        β = 50.0
        E_th = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        tol = max(1e-10, 10 * gap * exp(-β * gap))
        @test abs(E_th - E_gs) ≤ tol
    end
end

@testset "S1Heisenberg1D — direct ED cross-check at small N (β finite)" begin
    # Independent dense ED on the same matrix using the standard Boltzmann
    # average formula — guards against silent breakage in the eigen path.
    for (J, N, β) in ((1.0, 3, 1.0), (1.3, 4, 0.7), (0.6, 4, 2.5))
        m = S1Heisenberg1D(; J=J)
        H = QAtlas._s1_heisenberg_hamiltonian_matrix(m, N)
        evals = real.(eigvals(Hermitian(H)))
        emin = minimum(evals)
        ws = exp.(-β .* (evals .- emin))
        E_direct = sum(evals .* ws) / sum(ws)
        E_fetch = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        @test E_fetch ≈ E_direct rtol = 1e-12
    end
end

@testset "S1Heisenberg1D — ε = f + T·s identity (per-site)" begin
    # Self-validation across (Energy, FreeEnergy, ThermalEntropy).
    for (J, N, β) in ((1.0, 3, 0.5), (1.0, 4, 1.0), (1.3, 4, 2.0), (0.7, 5, 1.5))
        m = S1Heisenberg1D(; J=J)
        E_total = QAtlas.fetch(m, Energy(), OBC(N); beta=β)
        f = QAtlas.fetch(m, FreeEnergy(), OBC(N); beta=β)
        s = QAtlas.fetch(m, ThermalEntropy(), OBC(N); beta=β)
        ε = E_total / N
        @test isapprox(ε, f + s / β; atol=1e-10, rtol=1e-10)
    end
end

@testset "S1Heisenberg1D — c = β² Var(H)/N matches AutoDiff -β² ∂ε/∂β" begin
    # The variance formula for c is exact (no AD); cross-check it against
    # the AD path that PR #115 established for TFIM.
    using ForwardDiff
    for (J, N, β) in ((1.0, 3, 0.5), (1.0, 4, 1.0), (1.3, 4, 2.0))
        m = S1Heisenberg1D(; J=J)
        c_direct = QAtlas.fetch(m, SpecificHeat(), OBC(N); beta=β)
        dE_dβ = ForwardDiff.derivative(b -> QAtlas.fetch(m, Energy(), OBC(N); beta=b), β)
        c_ad = -β^2 * dE_dβ / N
        @test isapprox(c_direct, c_ad; atol=1e-9, rtol=1e-9)
    end
end

@testset "S1Heisenberg1D — monotone cooling E(β₁) ≥ E(β₂) for β₁ < β₂" begin
    m = S1Heisenberg1D(; J=1.0)
    Es = [QAtlas.fetch(m, Energy(), OBC(4); beta=β) for β in (0.0, 0.5, 1.0, 2.0, 5.0)]
    @test issorted(Es; rev=true)
    @test isapprox(Es[1], 0.0; atol=1e-12)
end

@testset "S1Heisenberg1D — dense ED cap" begin
    @test_throws ArgumentError QAtlas._spin1_string(
        QAtlas._MAX_ED_SITES_S1 + 1, 1 => QAtlas._S1_x
    )
    @test_throws ArgumentError QAtlas._s1_heisenberg_hamiltonian_matrix(S1Heisenberg1D(), 1)
end
