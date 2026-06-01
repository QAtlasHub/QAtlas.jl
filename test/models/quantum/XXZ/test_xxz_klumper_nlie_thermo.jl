# test/models/quantum/XXZ/test_xxz_klumper_nlie_thermo.jl
#
# Part of the Klümper-NLIE-backed FreeEnergy@Infinite test suite for
# `XXZ1D`, split from the original `test_xxz_klumper_nlie.jl` so each
# subset runs in its own CI shard.

using Test
using QAtlas
using QAtlas: XXZ1D, FreeEnergy, ThermalEntropy, SpecificHeat, Energy, Infinite

using LinearAlgebra: Hermitian, eigen
using SparseArrays: spzeros, sparse, kron
using QAtlas: _XXZ_NLIE_GRID_CACHE

@testset "XXZ Klümper NLIE — Δ = 0.5 thermodynamic family" begin
    # All Δ = 0.5 testsets share one cached grid (γ = π/3, built ~90 s on
    # first call). Group them together so the build is amortised across
    # all assertions in this file.

    @testset "high-T limit at Δ = 0.5: f − e₀ → −T·ln 2" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        β = 0.001
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        e0 = QAtlas.fetch(m, Energy{:per_site}(), Infinite())
        @test isfinite(f)
        @test isapprox(f - e0, -log(2) / β; rtol=1e-3)
    end

    @testset "Δ = 0.5 matches ED@N = 10 within boundary-CFT budget" begin
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
        @test isapprox(f_nlie, f_ed; rtol=0.05)
    end

    @testset "grid cache hit on repeated β" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0)
        n_before = length(_XXZ_NLIE_GRID_CACHE)
        QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=2.0)
        n_after = length(_XXZ_NLIE_GRID_CACHE)
        @test n_after == n_before
    end

    @testset "f(β) is monotone non-decreasing in β at Δ = 0.5" begin
        m = XXZ1D(; J=1.0, Δ=0.5)
        βs = (0.1, 0.5, 1.0, 2.0, 5.0)
        fs = [QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β) for β in βs]
        @test all(isfinite, fs)
        @test all(diff(fs) .>= -1e-6)
    end

    @testset "Δ = -0.5 (γ > π/2 half) — TODO: needs γ > π/2 solver fix" begin
        m = XXZ1D(; J=1.0, Δ=-0.5)
        @test_skip isfinite(QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=1.0))
    end
end
