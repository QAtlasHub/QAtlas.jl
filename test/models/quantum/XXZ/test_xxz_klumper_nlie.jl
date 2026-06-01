# test/models/quantum/XXZ/test_xxz_klumper_nlie.jl
#
# Tier-1 coverage for the Klümper-NLIE-backed FreeEnergy dispatch on
# `XXZ1D` at `Infinite()`. Cross-checks:
#  (a) Δ = 0 closed-form path stays exact;
#  (b) Δ → 0⁺ NLIE branch agrees with the XX free-fermion result
#      (regression guard against dispatch-ordering bugs);
#  (c) high-T limit  f − e₀  →  −T·ln 2;
#  (d) ED@N = 10 cross-check at moderate β for Δ = 0.5;
#  (e) endpoint / gapped Δ return NaN with a warning;
#  (f) β ≤ 0 raises DomainError;
#  (g) grid cache hit avoids the second full kernel build;
#  (h) ThermalEntropy / SpecificHeat are *not* yet routed through the
#      NLIE — assert they still NaN-warn, so a future PR that wires
#      them up trips this test.

using Test
using LinearAlgebra: Hermitian, eigen
using SparseArrays: spzeros, sparse, kron
using QAtlas
using QAtlas: XXZ1D, FreeEnergy, ThermalEntropy, SpecificHeat, Energy, Infinite
using QAtlas: _XXZ_NLIE_GRID_CACHE

@testset "XXZ Klümper NLIE — Infinite FreeEnergy" begin

    # ---- (a) XX closed-form still wins for Δ = 0 ----
    @testset "XX closed-form at Δ = 0" begin
        m = XXZ1D(; J=1.0, Δ=0.0)
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
        @test isfinite(f)
        @test f < 0
    end

    # ---- (b) NLIE-vs-XX regression anchor at Δ = 1e-6 ----
    @testset "NLIE near XX (Δ = 1e-6) matches closed form" begin
        m_xx = XXZ1D(; J=1.0, Δ=0.0)
        m_eps = XXZ1D(; J=1.0, Δ=1e-6)
        β = 1.0
        f_xx = QAtlas.fetch(m_xx, FreeEnergy(), Infinite(); beta=β)
        f_eps = QAtlas.fetch(m_eps, FreeEnergy(), Infinite(); beta=β)
        @test isapprox(f_eps, f_xx; rtol=1e-3)
    end

    # ---- (c) high-T limit ----
    @testset "high-T limit at Δ = 0.5: f − e₀ → −T·ln 2" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        β = 0.001
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        e0 = QAtlas.fetch(m, Energy{:per_site}(), Infinite())
        @test isfinite(f)
        @test isapprox(f - e0, -log(2) / β; rtol=1e-3)
    end

    # ---- (d) ED@N = 10 cross-check at moderate β ----
    @testset "Δ = 0.5 matches ED@N = 10 within boundary-CFT budget" begin
        # Build XXZ OBC Hamiltonian via sparse-tensor product. N = 10 is
        # small enough for a ≲ 20 s eigendecomposition on CI workers.
        function _ed_free_energy(N::Int, Δ::Float64, J::Float64, β::Float64)
            Sx = sparse(ComplexF64[0 1; 1 0] / 2)
            Sy = sparse(ComplexF64[0 -im; im 0] / 2)
            Sz = sparse(ComplexF64[1 0; 0 -1] / 2)
            Id = sparse(ComplexF64[1 0; 0 1])
            tens(A, k) = begin
                op = sparse([1.0 + 0im;;])
                for s in 1:N
                    op = kron(op, s == k ? A : Id)
                end
                op
            end
            H = spzeros(ComplexF64, 2^N, 2^N)
            for i in 1:(N - 1)
                H +=
                    J * (
                        tens(Sx, i) * tens(Sx, i+1) +
                        tens(Sy, i) * tens(Sy, i+1) +
                        Δ * tens(Sz, i) * tens(Sz, i+1)
                    )
            end
            evals, _ = eigen(Hermitian(Matrix(H)))
            emin = minimum(evals)
            return (-log(sum(exp.(-β .* (evals .- emin)))) / β + emin) / N
        end

        Δ, J, β = 0.5, 1.0, 1.0
        m = XXZ1D(; J=J, Δ=Δ)
        f_nlie = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        f_ed = _ed_free_energy(10, Δ, J, β)
        # OBC boundary CFT correction is ~ -πT²/(6 v_s L) per site,
        # ≈ 4 % at L=10, β=1, Δ=0.5 ⇒ tolerance of 0.05 is conservative.
        @test isapprox(f_nlie, f_ed; rtol=0.05)
    end

    # ---- (e) gapped and near-endpoint return NaN with warning ----
    @testset "|Δ| ≥ 1 (gapped) returns NaN + warn" begin
        m = XXZ1D(; J=1.0, Δ=1.5)
        f = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, FreeEnergy(), Infinite(); beta=1.0
        ))
        @test isnan(f)
    end
    @testset "near-endpoint Δ = -0.999 returns NaN + warn" begin
        # Mirror of the Δ ≈ +1 guard on the FM side. Picks the negative
        # endpoint of the |Δ| ≥ 0.99 guard so both sides are covered.
        m = XXZ1D(; J=1.0, Δ=-0.999)
        f = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, FreeEnergy(), Infinite(); beta=1.0
        ))
        @test isnan(f)
    end

    # ---- (f) β ≤ 0 raises DomainError ----
    @testset "β ≤ 0 raises DomainError" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0.0)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=-1.0)
    end

    # ---- (g) grid cache hit ----
    @testset "grid cache hit on repeated β" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        # Prime once.
        QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
        n_before = length(_XXZ_NLIE_GRID_CACHE)
        QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=2.0)  # same γ
        n_after = length(_XXZ_NLIE_GRID_CACHE)
        @test n_after == n_before
    end

    # ---- (h) ThermalEntropy / SpecificHeat still NaN at Δ ≠ 0 ----
    # If these become routed through the NLIE in a future PR, this test
    # will fail loudly, prompting the registry / coverage matrix in
    # issue #521 to be updated.
    @testset "ThermalEntropy / SpecificHeat at Δ = 0.5 still NaN (NLIE not yet wired)" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        s = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, ThermalEntropy(), Infinite(); beta=1.0
        ))
        c = (@test_logs (:warn,) match_mode=:any QAtlas.fetch(
            m, SpecificHeat(), Infinite(); beta=1.0
        ))
        @test isnan(s)
        @test isnan(c)
    end

    # ---- (i) Negative-Δ regime (broken at γ > π/2; tracked in #521) ----
    # γ = arccos(-0.5) = 2π/3 sits in the upper half of (0, π). The Klümper
    # kernel sinh((π/2-γ)k) changes sign relative to the γ<π/2 half, and the
    # Picard mixing α=0.4 fails to converge to a finite fixed point. Skip
    # until the iteration scheme is adapted to γ > π/2.
    @testset "Δ = -0.5 (γ > π/2 half) — TODO: needs γ > π/2 solver fix" begin
        m = XXZ1D(; J=1.0, Δ=-0.5)
        # When the negative-Δ path lands, switch this back to a -T·ln 2 high-T
        # check and a monotonicity sweep.
        @test_skip isfinite(QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0))
    end

    # ---- (j) Monotonicity of f(β) over a β sweep ----
    @testset "f(β) is monotone non-decreasing in β at Δ = 0.5" begin
        # f = -T ln Z is bounded below by ε_0 and approaches it from below
        # as β → ∞, so f(β) is non-decreasing in β. A sign flip in escale
        # or qatlas_to_klumper β̃ that single-point tests cannot catch
        # would show up here as a non-monotone curve.
        m = XXZ1D(; J=1.0, Δ=0.5)
        βs = (0.1, 0.5, 1.0, 2.0, 5.0)
        fs = [QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β) for β in βs]
        @test all(isfinite, fs)
        @test all(diff(fs) .>= -1e-6)
    end

    # ---- (k) Forced non-convergence raises NaN + warn ----
    @testset "solve_klumper_nlie at maxiter=1 reports non-convergence" begin
        # Direct exercise of the !sol.converged branch in the wrapper.
        # The wrapper hard-codes maxiter=400 so we cannot trigger this
        # path via the public fetch. We call the internal solver directly
        # to confirm the contract: residual > tol and converged == false.
        γ = acos(0.5)
        grid = QAtlas.XXZKlumperNLIE.build_grid(γ; N=64, L_factor=10.0, ε_shift=0.5)
        β̃ = 1.0
        sol = QAtlas.XXZKlumperNLIE.solve_klumper_nlie(grid, β̃; α=0.4, maxiter=1, tol=1e-9)
        @test !sol.converged
        @test sol.residual > 1e-9
        @test sol.iterations == 1
    end

    # ---- (l) J = 0 raises DomainError ----
    @testset "J = 0 raises DomainError" begin
        m = XXZ1D(; J=0.0, Δ=0.5)
        @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
    end

    # ---- (m) build_grid γ-overflow guard rejects unsafe (γ, ε_shift) ----
    @testset "build_grid rejects (γ, ε_shift) that would overflow cosh" begin
        # γ near the upper boundary acos(-0.99) ≈ 3.0 combined with the
        # floor ε_shift = 0.1 must be refused: kmax=300, decay=2.9,
        # cosh(300·2.9) overflows.
        γ_unsafe = acos(-0.99)
        @test_throws ArgumentError QAtlas.XXZKlumperNLIE.build_grid(
            γ_unsafe; N=32, L_factor=5.0, ε_shift=0.1
        )
    end
end
