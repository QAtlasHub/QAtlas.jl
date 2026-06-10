# =============================================================================
# Standalone tests for XXZ1D infinite-volume spinon kinematics and exact 
# two-spinon longitudinal dynamical structure factor S^{zz}(q, ω) (Pérez Castillo 2020).
# =============================================================================

using QAtlas, Test

@testset "XXZ1D exact two-spinon longitudinal DSF" begin
    # 1. Parameter validations
    @testset "Domain validations" begin
        # J <= 0 throws DomainError
        @test_throws DomainError QAtlas.fetch(XXZ1D(; J=-1.0, Δ=1.5), ZZStructureFactor(), Infinite(); q=1.0, ω=1.0)
        # Δ <= 1.0 throws DomainError
        @test_throws DomainError QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), ZZStructureFactor(), Infinite(); q=1.0, ω=1.0)
        @test_throws DomainError QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.5), ZZStructureFactor(), Infinite(); q=1.0, ω=1.0)
    end

    # 2. Support and kinematics
    @testset "Continuum support" begin
        model = XXZ1D(; J=1.0, Δ=1.5)
        # Choose a momentum q
        q = π / 2
        
        # Elliptic parameters to find the continuum boundaries
        qe, k, kp, K_val, Kp_val, ε = QAtlas._xxz_elliptic_params(model.Δ)
        κ = (1.0 - kp) / (1.0 + kp)
        I_val = (K_val / pi) * sqrt(model.Δ^2 - 1.0)
        
        omega0 = (2.0 * I_val / (1.0 + κ)) * sin(q)
        omegaminus = (2.0 * I_val / (1.0 + κ)) * sqrt(1.0 + κ^2 - 2.0 * κ * cos(q))
        
        # Test zero outside support
        @test QAtlas.fetch(model, ZZStructureFactor(), Infinite(); q=q, ω=0.5 * omega0) == 0.0
        @test QAtlas.fetch(model, ZZStructureFactor(), Infinite(); q=q, ω=1.5 * omegaminus) == 0.0
        
        # Test finite positive inside support
        ωmid = (omega0 + omegaminus) / 2
        S = QAtlas.fetch(model, ZZStructureFactor(), Infinite(); q=q, ω=ωmid)
        @test isfinite(S)
        @test S > 0
    end

    # 3. Limit checking to Heisenberg1D exact 2-spinon formula
    @testset "Isotropic limit Δ -> 1+" begin
        q = π / 2
        # For J=1.0, Heisenberg lower/upper edges at q=pi/2:
        # ε_L = pi/2 ≈ 1.570796
        # ε_U = pi * sin(pi/4) ≈ 2.22144
        # Choose ω in the middle of the continuum
        ω = 1.9
        
        S_heisenberg = QAtlas.fetch(Heisenberg1D(), ZZStructureFactor(), Infinite(); q=q, ω=ω, method=:exact_2spinon)
        
        # Compare with XXZ at Δ -> 1.0
        for Δ in (1.001, 1.0001, 1.00001)
            S_xxz = QAtlas.fetch(XXZ1D(; J=1.0, Δ=Δ), ZZStructureFactor(), Infinite(); q=q, ω=ω, method=:exact_2spinon)
            # Check convergence
            @test isapprox(S_xxz, S_heisenberg, rtol=1e-2)
        end
    end
end
