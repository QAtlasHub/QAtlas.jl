using QAtlas, Test, LinearAlgebra

# Heisenberg1D OBC delegators to XXZ1D(Δ = 1).  Every observable here
# should agree bit-for-bit with the corresponding XXZ1D(; J, Δ = 1.0)
# fetch — exercise the full surface so the delegator stays in sync
# with future XXZ1D additions.

@testset "Heisenberg1D OBC: scalar thermo matches XXZ1D(Δ=1)" begin
    Js = (1.0, 1.5)
    Ns = (4, 6)
    βs = (0.3, 1.0, 2.5)
    for J in Js, N in Ns, β in βs
        model_xxz = XXZ1D(; J=J, Δ=1.0)

        # Energy total
        E_heis = QAtlas.fetch(Heisenberg1D(), Energy(), OBC(N); beta=β, J=J)
        E_xxz = QAtlas.fetch(model_xxz, Energy(), OBC(N); beta=β)
        @test E_heis ≈ E_xxz rtol = 1e-12

        # FreeEnergy / ThermalEntropy / SpecificHeat (per site)
        for Q in (FreeEnergy(), ThermalEntropy(), SpecificHeat())
            v_heis = QAtlas.fetch(Heisenberg1D(), Q, OBC(N); beta=β, J=J)
            v_xxz = QAtlas.fetch(model_xxz, Q, OBC(N); beta=β)
            @test v_heis ≈ v_xxz rtol = 1e-12
        end
    end
end

@testset "Heisenberg1D OBC: magnetisations + susceptibilities match XXZ1D(Δ=1)" begin
    J, N, β = 1.5, 5, 0.7
    model_xxz = XXZ1D(; J=J, Δ=1.0)
    for Q in (
        MagnetizationX(),
        MagnetizationY(),
        MagnetizationZ(),
        SusceptibilityXX(),
        SusceptibilityYY(),
        SusceptibilityZZ(),
    )
        v_heis = QAtlas.fetch(Heisenberg1D(), Q, OBC(N); beta=β, J=J)
        v_xxz = QAtlas.fetch(model_xxz, Q, OBC(N); beta=β)
        @test v_heis ≈ v_xxz atol = 1e-12
    end
end

@testset "Heisenberg1D OBC: local observables match XXZ1D(Δ=1)" begin
    J, N, β = 1.0, 5, 1.0
    model_xxz = XXZ1D(; J=J, Δ=1.0)
    for Q in
        (MagnetizationXLocal(), MagnetizationYLocal(), MagnetizationZLocal(), EnergyLocal())
        v_heis = QAtlas.fetch(Heisenberg1D(), Q, OBC(N); beta=β, J=J)
        v_xxz = QAtlas.fetch(model_xxz, Q, OBC(N); beta=β)
        @test all(isapprox.(v_heis, v_xxz; atol=1e-12))
    end
end

@testset "Heisenberg1D OBC: two-point correlators match XXZ1D(Δ=1)" begin
    J, N, β = 1.0, 5, 1.5
    model_xxz = XXZ1D(; J=J, Δ=1.0)
    for CorrTy in (XXCorrelation, YYCorrelation, ZZCorrelation),
        mode in (:static, :connected)

        Q = CorrTy(; mode=mode)
        for i in 1:N, j in i:N
            v_heis = QAtlas.fetch(Heisenberg1D(), Q, OBC(N); beta=β, i=i, j=j, J=J)
            v_xxz = QAtlas.fetch(model_xxz, Q, OBC(N); beta=β, i=i, j=j)
            @test v_heis ≈ v_xxz atol = 1e-12
        end
    end
end

@testset "Heisenberg1D OBC: VonNeumannEntropy + RenyiEntropy match XXZ1D(Δ=1)" begin
    J, N, β = 1.0, 5, Inf
    model_xxz = XXZ1D(; J=J, Δ=1.0)
    for ℓ in (1, 2, 3, 4)
        v_heis = QAtlas.fetch(Heisenberg1D(), VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=β, J=J)
        v_xxz = QAtlas.fetch(model_xxz, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=β)
        @test v_heis ≈ v_xxz atol = 1e-10

        for α in (0.5, 2.0, 3.0)
            q = RenyiEntropy(α)
            v_heis = QAtlas.fetch(Heisenberg1D(), q, OBC(N); ℓ=ℓ, beta=β, J=J)
            v_xxz = QAtlas.fetch(model_xxz, q, OBC(N); ℓ=ℓ, beta=β)
            @test v_heis ≈ v_xxz atol = 1e-10
        end
    end
end

@testset "Heisenberg1D: MassGap matches XXZ1D(Δ=1)" begin
    for J in (1.0, 1.5), N in (4, 5)
        v_heis = QAtlas.fetch(Heisenberg1D(), MassGap(), OBC(N); J=J)
        v_xxz = QAtlas.fetch(XXZ1D(; J=J, Δ=1.0), MassGap(), OBC(N))
        @test v_heis ≈ v_xxz rtol = 1e-12
    end
    # Infinite (gapless Luttinger): both return 0
    v_heis_inf = QAtlas.fetch(Heisenberg1D(), MassGap(), Infinite())
    @test v_heis_inf == 0.0
end
