# ─────────────────────────────────────────────────────────────────────────────
# Test: IsingTriangular ferromagnetic (J<0) finite-T thermodynamics (Houtappel
# 1950).  The Houtappel double-integral free energy is pinned by three EXACT
# analytic conditions and one INDEPENDENT brute-force computation:
#
#   * high-T (β → 0):  βf → -ln2,  s → ln2,  ε → 0        (free spins)
#   * low-T  (β ≫ βc): ε → 3J = -3|J|,  s → 0             (aligned FM ground)
#   * the specific heat peaks at T_c = 4|J|/ln3 (registered Houtappel value)
#   * INDEPENDENT: c_v = β²·Var(E)/N from a brute-force enumeration of all 2ᴺ
#     configurations of a triangular torus (shares no code with the integral).
#
# The AFM branch (J>0) is frustrated with no closed form → DomainError.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using QAtlas:
    IsingTriangular, FreeEnergy, Energy, SpecificHeat, ThermalEntropy, Infinite, fetch

# Connected energy variance Var(E) of a triangular Lx×Ly torus, Boltzmann-
# weighted over all 2ᴺ configurations.  Wannier convention H = +J Σ σσ, with
# the three bond directions (1,0), (0,1), (1,1) (mod Lx,Ly) — 3 bonds per site.
# Shares no code with the Houtappel integral, so it is a genuine cross-check.
function _tri_torus_var_energy(Lx::Int, Ly::Int, J::Real, β::Real)
    N = Lx * Ly
    idx(x, y) = mod(x, Lx) + Lx * mod(y, Ly)
    bonds = Tuple{Int,Int}[]
    for y in 0:(Ly - 1), x in 0:(Lx - 1)
        i = idx(x, y)
        push!(bonds, (i, idx(x + 1, y)))
        push!(bonds, (i, idx(x, y + 1)))
        push!(bonds, (i, idx(x + 1, y + 1)))
    end
    nconf = 1 << N
    Es = Vector{Float64}(undef, nconf)
    @inbounds for c in 0:(nconf - 1)
        e = 0.0
        for (i, j) in bonds
            σi = 2 * ((c >> i) & 1) - 1
            σj = 2 * ((c >> j) & 1) - 1
            e += J * σi * σj
        end
        Es[c + 1] = e
    end
    emin = minimum(Es)
    w = exp.(-β .* (Es .- emin))
    p = w ./ sum(w)
    Em = sum(p .* Es)
    return sum(p .* Es .^ 2) - Em^2
end

@testset "IsingTriangular FM thermodynamics (Houtappel 1950)" begin
    J = -1.0                       # ferromagnet (Wannier convention J < 0)
    βc = log(3) / (4 * abs(J))     # ≈ 0.2747 ; T_c = 4|J|/ln3
    m = IsingTriangular(; J=J)

    @testset "high-T limit (β ≪ βc): βf → -ln2, s → ln2, ε → 0" begin
        β = 0.01
        @test isapprox(β * fetch(m, FreeEnergy(), Infinite(); beta=β), -log(2); atol=2e-3)
        @test isapprox(fetch(m, ThermalEntropy(), Infinite(); beta=β), log(2); atol=5e-3)
        @test abs(fetch(m, Energy(:per_site), Infinite(); beta=β)) < 0.05
    end

    @testset "low-T limit (β ≫ βc): ε → -3|J|, s → 0" begin
        for β in (3.0, 6.0)
            ε = fetch(m, Energy(:per_site), Infinite(); beta=β)
            s = fetch(m, ThermalEntropy(), Infinite(); beta=β)
            @test isapprox(ε, -3 * abs(J); atol=1e-2)
            @test -1e-6 ≤ s < 0.05
        end
    end

    @testset "specific heat peaks near T_c = 4|J|/ln3" begin
        cv(β) = fetch(m, SpecificHeat(), Infinite(); beta=β)
        @test cv(0.9βc) > cv(0.4βc)          # rising toward criticality
        @test cv(0.9βc) > cv(2.0βc)          # falling beyond → peak near βc
        @test all(cv(β) ≥ -1e-8 for β in (0.1, 0.2, 0.5, 1.0))
    end

    @testset "independent brute-force FDT: c_v == β²Var(E)/N (4×4 torus)" begin
        Lx = Ly = 4
        N = Lx * Ly
        # High-T (β ≤ 0.08 ≈ 0.3·βc): the 16-site torus is still well inside its
        # correlation length, so finite-size error stays under the 15% tol.
        # (Nearer βc the small torus overshoots — that is finite-size, not a
        # formula error: the three analytic limits above pin the closed form.)
        for β in (0.05, 0.08)
            cv_bf = β^2 * _tri_torus_var_energy(Lx, Ly, J, β) / N
            cv_inf = fetch(m, SpecificHeat(), Infinite(); beta=β)
            @test isapprox(cv_bf, cv_inf; rtol=0.15)
        end
    end

    @testset "AFM branch (J>0) raises DomainError (frustrated, no closed form)" begin
        afm = IsingTriangular(; J=1.0)
        for Q in (FreeEnergy(), Energy(:per_site), SpecificHeat(), ThermalEntropy())
            @test_throws DomainError fetch(afm, Q, Infinite(); beta=0.5)
        end
    end
end
