using QAtlas, Test, LinearAlgebra

# Spin-1 Heisenberg observables (PR Tier 1) — small-N dense-ED reference.
# Cross-checks every fetch method against an independent direct-ED
# evaluation of `Tr(O ρ_β)` on the 3^N Hilbert space.

# Bring private helpers into scope for direct-ED cross checks.
using QAtlas: _S1_x, _S1_y, _S1_z, _spin1_string, _s1_heisenberg_hamiltonian_matrix
using QAtlas: _s1_thermal_kernel, _s1_partial_trace_A

# ──────────────────────────────────────────────────────────────────────
# Helper: independent thermal density matrix (3^N × 3^N).
# ──────────────────────────────────────────────────────────────────────
function _direct_thermal_rho(H::AbstractMatrix, β::Real)
    F = eigen(Hermitian(H))
    if isinf(β) && β > 0
        emin = F.values[1]
        mask = (F.values .- emin) .≤ 1e-12
        ws = zeros(Float64, length(F.values))
        ws[mask] .= 1.0 / count(mask)
    else
        emin = minimum(F.values)
        ws = exp.(-β .* (F.values .- emin))
        ws ./= sum(ws)
    end
    return F.vectors * Diagonal(ws) * F.vectors'
end
_direct_expect(O, ρ) = real(tr(O * ρ))

# ──────────────────────────────────────────────────────────────────────
# T → ∞ (β → 0) limits
# ──────────────────────────────────────────────────────────────────────

@testset "S1Heisenberg1D — β → 0:  ⟨Sᵅ⟩ → 0  (each Sᵅ traceless)" begin
    m = S1Heisenberg1D(; J=1.0)
    for N in (3, 4)
        β = 1e-4  # small but nonzero so ForwardDiff is happy
        @test abs(QAtlas.fetch(m, MagnetizationX(), OBC(N); beta=β)) < 1e-10
        @test abs(QAtlas.fetch(m, MagnetizationY(), OBC(N); beta=β)) < 1e-10
        @test abs(QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β)) < 1e-10
    end
end

@testset "S1Heisenberg1D — β → 0:  χ_αα → 2β/3 per site" begin
    # `Tr((Sᵅ)²) / 3 = 2/3` for the 3-dim Sᵅ matrices (spin-1, eigenvalues ±1, 0).
    # Across N independent traceless sites at β → 0, ⟨M_α⟩ → 0 and
    # ⟨M_α²⟩ → N · 2/3, so χ_αα = β · (2N/3) / N = 2β/3.
    m = S1Heisenberg1D(; J=1.0)
    β = 1e-4
    for N in (3, 4)
        for Q in (SusceptibilityXX(), SusceptibilityYY(), SusceptibilityZZ())
            χ = QAtlas.fetch(m, Q, OBC(N); beta=β)
            @test isapprox(χ, 2β / 3; atol=1e-8, rtol=1e-6)
        end
    end
end

# ──────────────────────────────────────────────────────────────────────
# T → 0 (β large): each thermal observable converges to the GS expectation.
# ──────────────────────────────────────────────────────────────────────

@testset "S1Heisenberg1D — β = 50 thermal vs GS expectation" begin
    for (J, N) in ((1.0, 3), (0.7, 4))
        m = S1Heisenberg1D(; J=J)
        H = _s1_heisenberg_hamiltonian_matrix(m, N)
        F = eigen(Hermitian(H))
        # Use the same equal-weight GS-projector convention the kernel uses
        # at β = Inf (avoids degeneracy mismatch with a single eigenvector).
        ρ_gs = _direct_thermal_rho(H, Inf)
        β = 50.0
        gap = F.values[2] - F.values[1]
        # tolerance: residue ≲ exp(-β · gap), keep one safety decade.
        tol = max(1e-9, 100 * exp(-β * gap))

        for (Q, O) in (
            (MagnetizationX(), sum(_spin1_string(N, i => _S1_x) for i in 1:N) / N),
            (MagnetizationY(), sum(_spin1_string(N, i => _S1_y) for i in 1:N) / N),
            (MagnetizationZ(), sum(_spin1_string(N, i => _S1_z) for i in 1:N) / N),
        )
            v_th = QAtlas.fetch(m, Q, OBC(N); beta=β)
            v_gs = _direct_expect(O, ρ_gs)
            @test isapprox(v_th, v_gs; atol=tol)
        end
    end
end

# ──────────────────────────────────────────────────────────────────────
# Direct-ED parity at small N for every fetch (susceptibility, correlators,
# local quantities, MassGap, entanglement)
# ──────────────────────────────────────────────────────────────────────

@testset "S1Heisenberg1D — N=2 dimer + finite β: direct ED matches fetch" begin
    # Dimer (N=2) — 9-state spectrum is closed-form ([-2J, -J×3, +J×5]).
    # Build everything from scratch to catch silent regressions.
    for (J, β) in ((1.0, 0.5), (1.3, 1.0), (0.7, 2.0))
        m = S1Heisenberg1D(; J=J)
        N = 2
        H = _s1_heisenberg_hamiltonian_matrix(m, N)
        ρ = _direct_thermal_rho(H, β)
        Mx = sum(_spin1_string(N, i => _S1_x) for i in 1:N)
        Mz = sum(_spin1_string(N, i => _S1_z) for i in 1:N)
        # Variance susceptibility direct
        χxx_direct = β * (_direct_expect(Mx * Mx, ρ) - _direct_expect(Mx, ρ)^2) / N
        χzz_direct = β * (_direct_expect(Mz * Mz, ρ) - _direct_expect(Mz, ρ)^2) / N
        @test QAtlas.fetch(m, SusceptibilityXX(), OBC(N); beta=β) ≈ χxx_direct rtol=1e-12
        @test QAtlas.fetch(m, SusceptibilityZZ(), OBC(N); beta=β) ≈ χzz_direct rtol=1e-12

        # Static + connected correlators at (i, j) = (1, 2)
        Sz1 = _spin1_string(N, 1 => _S1_z)
        Sz2 = _spin1_string(N, 2 => _S1_z)
        zz_static = _direct_expect(Sz1 * Sz2, ρ)
        zz_connected = zz_static - _direct_expect(Sz1, ρ) * _direct_expect(Sz2, ρ)
        @test QAtlas.fetch(m, ZZCorrelation(; mode=:static), OBC(N); beta=β, i=1, j=2) ≈
            zz_static rtol=1e-12
        @test QAtlas.fetch(m, ZZCorrelation(; mode=:connected), OBC(N); beta=β, i=1, j=2) ≈
            zz_connected rtol=1e-12

        # i = j case: ⟨Sz²⟩ on a single site
        Sz2_local = Sz1 * Sz1
        zz_ii = _direct_expect(Sz2_local, ρ)
        @test QAtlas.fetch(m, ZZCorrelation(; mode=:static), OBC(N); beta=β, i=1, j=1) ≈
            zz_ii rtol=1e-12
    end
end

@testset "S1Heisenberg1D — XX/YY correlators at small N" begin
    for (J, N, β) in ((1.0, 3, 0.6), (1.3, 4, 1.0))
        m = S1Heisenberg1D(; J=J)
        H = _s1_heisenberg_hamiltonian_matrix(m, N)
        ρ = _direct_thermal_rho(H, β)
        for (Q, S) in
            ((XXCorrelation(; mode=:static), _S1_x), (YYCorrelation(; mode=:static), _S1_y))
            for (i, j) in ((1, N), (1, 2), (2, 2))
                O = i == j ? _spin1_string(N, i => S * S) : _spin1_string(N, i => S, j => S)
                expected = _direct_expect(O, ρ)
                @test QAtlas.fetch(m, Q, OBC(N); beta=β, i=i, j=j) ≈ expected rtol=1e-10
            end
        end
    end
end

# ──────────────────────────────────────────────────────────────────────
# MassGap, OBC + Infinite
# ──────────────────────────────────────────────────────────────────────

@testset "S1Heisenberg1D — MassGap(OBC) = E₁ - E₀" begin
    for (J, N) in ((1.0, 4), (0.7, 5), (1.3, 3))
        m = S1Heisenberg1D(; J=J)
        H = _s1_heisenberg_hamiltonian_matrix(m, N)
        evals = sort(eigvals(Hermitian(H)))
        @test QAtlas.fetch(m, MassGap(), OBC(N)) ≈ evals[2] - evals[1] rtol=1e-12
    end
end

@testset "S1Heisenberg1D — Infinite-system literature constants" begin
    # Per-site GS energy density — White-Huse 1993 DMRG value.
    @test QAtlas.fetch(S1Heisenberg1D(; J=1.0), Energy(:per_site), Infinite()) ≈
        -1.40148403897 rtol=1e-12
    @test QAtlas.fetch(S1Heisenberg1D(; J=1.7), Energy(:per_site), Infinite()) ≈
        -1.40148403897 * 1.7 rtol=1e-12

    # Haldane gap — White-Huse 1993 DMRG value.
    @test QAtlas.fetch(S1Heisenberg1D(; J=1.0), MassGap(), Infinite()) ≈ 0.41048 rtol=1e-12
    @test QAtlas.fetch(S1Heisenberg1D(; J=2.5), MassGap(), Infinite()) ≈ 0.41048 * 2.5 rtol=1e-12
end

# ──────────────────────────────────────────────────────────────────────
# Local observables — sum identity
# ──────────────────────────────────────────────────────────────────────

@testset "S1Heisenberg1D — local observables sum to bulk values" begin
    for (J, N, β) in ((1.0, 3, 0.5), (1.3, 4, 1.0), (0.7, 4, 2.0))
        m = S1Heisenberg1D(; J=J)
        # Mz_local sum = N * MagnetizationZ
        mz_local = QAtlas.fetch(m, MagnetizationZLocal(), OBC(N); beta=β)
        mz_bulk = QAtlas.fetch(m, MagnetizationZ(), OBC(N); beta=β)
        @test sum(mz_local) ≈ N * mz_bulk atol=1e-10
        # Same for X
        mx_local = QAtlas.fetch(m, MagnetizationXLocal(), OBC(N); beta=β)
        mx_bulk = QAtlas.fetch(m, MagnetizationX(), OBC(N); beta=β)
        @test sum(mx_local) ≈ N * mx_bulk atol=1e-10
        # Energy_local sum = ⟨H⟩
        ε_local = QAtlas.fetch(m, EnergyLocal(), OBC(N); beta=β)
        E_total = QAtlas.fetch(m, Energy(:total), OBC(N); beta=β)
        @test sum(ε_local) ≈ E_total atol=1e-10
        @test length(ε_local) == N
        @test length(mx_local) == N
        @test length(mz_local) == N
    end
end

# ──────────────────────────────────────────────────────────────────────
# Entanglement entropies
# ──────────────────────────────────────────────────────────────────────

@testset "S1Heisenberg1D — VonNeumannEntropy basics (ground state)" begin
    # Haldane chain at OBC has quasi-degenerate GS manifold from end-state
    # physics (`_s1_thermal_kernel(model, N, Inf)` mixes the manifold with
    # equal weights once the splitting falls below 1e-12).  The resulting ρ
    # is a mixed state, so the pure-state Schmidt identity `S(ρ_A) =
    # S(ρ_{A^c})` does *not* apply.  We check only the universally valid
    # bounds:
    #   * S ≥ 0
    #   * S(ρ_A) ≤ log(dim A) = ℓ · log 3  (max entropy on the subsystem)
    for (J, N) in ((1.0, 4), (1.0, 5))
        m = S1Heisenberg1D(; J=J)
        Ss = [QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ) for ℓ in 1:(N - 1)]
        @test all(S -> S ≥ -1e-12, Ss)
        for ℓ in 1:(N - 1)
            @test Ss[ℓ] ≤ ℓ * log(3) + 1e-10
        end
    end
end

@testset "S1Heisenberg1D — VonNeumannEntropy at β → 0 saturates ℓ log 3" begin
    # Maximally mixed thermal state ρ ∝ I gives ρ_A = I/3^ℓ → S = ℓ log 3.
    m = S1Heisenberg1D(; J=1.0)
    for N in (3, 4), ℓ in 1:(N - 1)
        S = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=1e-6)
        @test isapprox(S, ℓ * log(3); atol=1e-3)
    end
end

@testset "S1Heisenberg1D — RenyiEntropy reduces to VonNeumann at α → 1" begin
    # log Tr ρ^α / (1-α) → -Tr ρ log ρ as α → 1.
    m = S1Heisenberg1D(; J=1.0)
    N = 4
    ℓ = 2
    S_vN = QAtlas.fetch(m, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=1.0)
    for α in (0.99, 1.01)
        S_α = QAtlas.fetch(m, RenyiEntropy(α), OBC(N); ℓ=ℓ, beta=1.0)
        @test isapprox(S_α, S_vN; atol=5e-2)
    end
    # Strict inequality S_2 ≤ S_vN for α > 1 (Rényi is non-increasing in α)
    S_2 = QAtlas.fetch(m, RenyiEntropy(2), OBC(N); ℓ=ℓ, beta=1.0)
    @test S_2 ≤ S_vN + 1e-12
    # And S_{1/2} ≥ S_vN for α < 1
    S_half = QAtlas.fetch(m, RenyiEntropy(0.5), OBC(N); ℓ=ℓ, beta=1.0)
    @test S_half ≥ S_vN - 1e-12
end

# ──────────────────────────────────────────────────────────────────────
# Energy{:per_site} <-> Energy{:total} routing for OBC native :total
# ──────────────────────────────────────────────────────────────────────

@testset "S1Heisenberg1D — Energy(:per_site) routes through :total / N at OBC" begin
    m = S1Heisenberg1D(; J=1.0)
    for (N, β) in ((3, 0.5), (4, 1.0))
        E_total = QAtlas.fetch(m, Energy(:total), OBC(N); beta=β)
        ε = QAtlas.fetch(m, Energy(:per_site), OBC(N); beta=β)
        @test ε ≈ E_total / N rtol=1e-12
    end
end
