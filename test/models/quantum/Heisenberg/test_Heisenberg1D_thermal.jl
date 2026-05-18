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

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "Heisenberg1D thermal OBC — verification cards" begin
    Sx, Sy, Sz = spin_ops(1 // 2)

    # MassGap in the gapless Luttinger liquid regime is exactly 0
    verify(
        Heisenberg1D(),
        MassGap(),
        Infinite();
        route=:second_closed_form,
        independent=0.0,
        agree_within=1e-14,
        refs=["Heisenberg chain is gapless (des Cloizeaux-Pearson 1962): gap = 0"],
    )

    # OBC thermal energy at N=4, beta=1 vs independent ED (thermo_from_spectrum)
    let J = 1.0, N = 4, beta = 1.0
        bond = J * (kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz))
        E_ind, _, _, _ = thermo_from_spectrum(dense_spectrum(chain_hamiltonian(2, N, bond)), beta)
        verify(
            Heisenberg1D(),
            Energy(),
            OBC(N);
            route=:ed_finite_size,
            fetch_kw=(; beta=beta, J=J),
            independent=E_ind,
            agree_within=1e-9,
            refs=["Direct OBC ED via generic_ed chain_hamiltonian + thermo_from_spectrum"],
        )
    end

    # Delegation invariant: Heisenberg1D thermal OBC === XXZ1D(Delta=1) at same J
    let beta = 1.0, N = 4
        verify(
            Heisenberg1D(),
            Energy(),
            OBC(N);
            route=:delegation_invariant,
            fetch_kw=(; beta=beta, J=1.5),
            independent=QAtlas.fetch(XXZ1D(; J=1.5, Δ=1.0), Energy(), OBC(N); beta=beta),
            agree_within=1e-12,
            refs=["Heisenberg1D thermal OBC delegates to XXZ1D(Delta=1): same J must match"],
        )
    end
end
