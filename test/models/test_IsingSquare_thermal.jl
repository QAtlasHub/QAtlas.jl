using QAtlas, Test, LinearAlgebra
using ForwardDiff: ForwardDiff

# ─── ED reference helpers ─────────────────────────────────────────────
# Brute-force the 2^(Lx·Ly) configuration sum for tiny lattices to
# build an exact `log Z` reference; ε / s / c_v then follow by the
# usual Boltzmann moments.

function _bruteforce_thermo(Lx::Int, Ly::Int, β::Real, J::Real=1.0)
    N = Lx * Ly
    energies = Float64[]
    for cfg in 0:(2^N - 1)
        σ = [(cfg >> k) & 1 == 1 ? 1 : -1 for k in 0:(N - 1)]
        idx(i, j) = (i - 1) * Ly + j  # row i, col j
        E = 0.0
        @inbounds for i in 1:Lx, j in 1:Ly
            ip = i % Lx + 1
            jp = j % Ly + 1
            E += -J * σ[idx(i, j)] * σ[idx(ip, j)]   # vertical (transfer dir)
            E += -J * σ[idx(i, j)] * σ[idx(i, jp)]   # horizontal
        end
        push!(energies, E)
    end
    e_min = minimum(energies)
    ws = exp.(-β .* (energies .- e_min))
    Z = sum(ws)
    log_Z = -β * e_min + log(Z)
    ε_per = sum(energies .* ws) / Z / N
    f_per = -log_Z / (β * N)
    s_per = β * (ε_per - f_per)
    H2 = sum(energies .^ 2 .* ws) / Z
    H1 = sum(energies .* ws) / Z
    c_per = β^2 * (H2 - H1^2) / N
    return (; ε_per, f_per, s_per, c_per, log_Z)
end

@testset "IsingSquare PBC: ED comparison at small (Lx, Ly)" begin
    J = 1.0
    for (Lx, Ly) in [(3, 3), (3, 4), (4, 3)]   # avoid Lx,Ly=2 (transfer-matrix
                                              # PBC double-counts bonds there)
        for β in (0.1, 0.3, 0.6)
            m = IsingSquare(; J=J)
            ed = _bruteforce_thermo(Lx, Ly, β, J)

            f = QAtlas.fetch(m, FreeEnergy(), PBC(0); beta=β, Lx=Lx, Ly=Ly)
            ε = QAtlas.fetch(m, Energy(:per_site), PBC(0); beta=β, Lx=Lx, Ly=Ly)
            s = QAtlas.fetch(m, ThermalEntropy(), PBC(0); beta=β, Lx=Lx, Ly=Ly)
            c = QAtlas.fetch(m, SpecificHeat(), PBC(0); beta=β, Lx=Lx, Ly=Ly)

            @test f ≈ ed.f_per atol = 1e-10
            @test ε ≈ ed.ε_per atol = 1e-7
            @test s ≈ ed.s_per atol = 1e-7
            @test c ≈ ed.c_per atol = 1e-3   # 2nd derivative central diff
        end
    end
end

@testset "IsingSquare PBC: Gibbs identity ε = f + T·s" begin
    m = IsingSquare(; J=1.0)
    for β in (0.1, 0.3, 0.6), (Lx, Ly) in [(3, 3), (4, 4)]
        f = QAtlas.fetch(m, FreeEnergy(), PBC(0); beta=β, Lx=Lx, Ly=Ly)
        ε = QAtlas.fetch(m, Energy(:per_site), PBC(0); beta=β, Lx=Lx, Ly=Ly)
        s = QAtlas.fetch(m, ThermalEntropy(), PBC(0); beta=β, Lx=Lx, Ly=Ly)
        @test ε ≈ f + s / β atol = 1e-10
    end
end

@testset "IsingSquare Infinite: high-T and low-T limits" begin
    m = IsingSquare(; J=1.0)

    # β → 0 (T → ∞): f → -log(2)/β, ε → 0, s → log(2)
    β_hi = 1e-4
    f_hi = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β_hi)
    ε_hi = QAtlas.fetch(m, Energy(:per_site), Infinite(); beta=β_hi)
    s_hi = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β_hi)
    @test f_hi ≈ -log(2) / β_hi atol = 1e-2 / β_hi    # tolerance scales with f
    @test abs(ε_hi) < 1e-3
    @test s_hi ≈ log(2) atol = 1e-3

    # β → ∞ (T → 0): per-site GS energy = -2J (every bond aligned, 2 bonds/site).
    β_lo = 50.0
    ε_lo = QAtlas.fetch(m, Energy(:per_site), Infinite(); beta=β_lo)
    @test ε_lo ≈ -2.0 atol = 1e-6
end

@testset "IsingSquare Infinite: PBC → Infinite at large Lx=Ly" begin
    m = IsingSquare(; J=1.0)
    β = 0.3   # below T_c (β_c ≈ 0.4407)
    f_inf = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
    f_pbc = QAtlas.fetch(m, FreeEnergy(), PBC(0); beta=β, Lx=8, Ly=8)
    @test abs(f_pbc - f_inf) < 1e-2
end

@testset "IsingSquare Infinite: Onsager c_v finite away from T_c" begin
    m = IsingSquare(; J=1.0)
    # away from T_c (β_c ≈ 0.4407)
    for β in (0.1, 0.3, 0.6, 1.0, 2.0)
        c = QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=β)
        @test isfinite(c)
        @test c > 0
    end
end

@testset "IsingSquare Infinite: SpecificHeat consistent with ForwardDiff on E" begin
    # The fetch implementation already uses central differences; this
    # cross-checks the closed form against an independent ForwardDiff
    # path on the energy.
    m = IsingSquare(; J=1.0)
    for β in (0.2, 0.6)
        c_fetch = QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=β)
        # c_v = -β² ∂ε/∂β
        dε = ForwardDiff.derivative(
            b -> QAtlas.fetch(m, Energy(:per_site), Infinite(); beta=b), β
        )
        c_ad = -β^2 * dε
        @test c_fetch ≈ c_ad rtol = 1e-2  # central-diff vs ForwardDiff slack
    end
end
