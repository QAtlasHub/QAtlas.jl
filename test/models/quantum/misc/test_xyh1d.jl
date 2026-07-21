# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: XYh1D — general anisotropic XY model exact solutions
# (Lieb-Schultz-Mattis 1961; Pfeuty 1970).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

if !isdefined(Main, :verify)
    include("../../../util/verify.jl")
end

@testset "XYh1D — MassGap (Infinite and OBC)" begin
    # h = 0: gapless XX chain
    @test QAtlas.fetch(XYh1D(), MassGap(), Infinite()) == 0.0
    # Inside band (|h| < 2J): gapless
    @test QAtlas.fetch(XYh1D(; h=1.0), MassGap(), Infinite()) == 0.0
    # At the critical h = 2J: gap closes
    @test QAtlas.fetch(XYh1D(; h=2.0), MassGap(), Infinite()) == 0.0
    # Polarised (|h| > 2J): finite gap = 2(|h| - 2J)
    @test QAtlas.fetch(XYh1D(; h=3.0), MassGap(), Infinite()) == 2.0
    @test QAtlas.fetch(XYh1D(; h=-3.0), MassGap(), Infinite()) == 2.0
    @test QAtlas.fetch(XYh1D(; h=5.0), MassGap(), Infinite()) == 6.0

    # Anisotropic Jx ≠ Jy (MassGap should be finite for h < Jx+Jy if Jx ≠ Jy)
    # e.g., Jx=1.0, Jy=0.5, h=0.0: MassGap = 2 * |Jx - Jy| = 2.0 * 0.5 = 1.0
    @test QAtlas.fetch(XYh1D(; Jx=1.0, Jy=0.5, h=0.0), MassGap(), Infinite()) ≈ 1.0

    # OBC MassGap matches the smallest positive eigenvalue of the BdG spectrum
    let m = XYh1D(; Jx=1.0, Jy=0.5, h=0.5)
        gap_inf = QAtlas.fetch(m, MassGap(), Infinite())
        gap_obc = QAtlas.fetch(m, MassGap(), OBC(200))
        @test isapprox(gap_obc, gap_inf; atol=1e-3)
    end
end

@testset "XYh1D — rejects Jx, Jy ≤ 0" begin
    @test_throws DomainError XYh1D(; Jx=0.0)
    @test_throws DomainError XYh1D(; Jx=-1.0)
    @test_throws DomainError XYh1D(; Jy=0.0)
    @test_throws DomainError XYh1D(; Jy=-1.0)
end

@testset "XYh1D — Energy{:per_site} (Infinite and OBC)" begin
    # h = 0, J = 1: E/N = -4/π  (Lieb-Schultz-Mattis 1961)
    e0 = QAtlas.fetch(XYh1D(), Energy{:per_site}(), Infinite())
    @test isapprox(e0, -4 / π; atol=1e-12)

    # J linearity at h = 0
    for J in (0.5, 1.0, 2.0, 3.5)
        e = QAtlas.fetch(XYh1D(; Jx=J, Jy=J, h=0.0), Energy{:per_site}(), Infinite())
        @test isapprox(e, -4J / π; atol=1e-12)
    end
end

@testset "XYh1D — Finite Temperature potentials" begin
    # Check that OBC free energy approaches the infinite thermodynamic limit
    let m = XYh1D(; Jx=1.2, Jy=0.8, h=0.6), β = 2.0
        f_inf = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        f_obc = QAtlas.fetch(m, FreeEnergy(), OBC(200); beta=β)
        @test isapprox(f_obc, f_inf; atol=1e-2)

        s_inf = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β)
        s_obc = QAtlas.fetch(m, ThermalEntropy(), OBC(200); beta=β)
        @test isapprox(s_obc, s_inf; atol=1e-2)

        c_inf = QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=β)
        c_obc = QAtlas.fetch(m, SpecificHeat(), OBC(200); beta=β)
        @test isapprox(c_obc, c_inf; atol=1e-2)
    end
end

@testset "XYh1D — MagnetizationZ and SusceptibilityZZ" begin
    let m = XYh1D(; Jx=1.2, Jy=0.8, h=0.6), β = 2.0
        mz_inf = QAtlas.fetch(m, MagnetizationZ(), Infinite(); beta=β)
        mz_obc = QAtlas.fetch(m, MagnetizationZ(), OBC(200); beta=β)
        @test isapprox(mz_obc, mz_inf; atol=1e-2)

        χ_inf = QAtlas.fetch(m, SusceptibilityZZ(), Infinite(); beta=β)
        χ_obc = QAtlas.fetch(m, SusceptibilityZZ(), OBC(200); beta=β)
        @test isfinite(χ_inf) && χ_inf > 0
        @test isfinite(χ_obc) && χ_obc > 0
        # NOTE: OBC and Infinite normalisations differ by a factor ~2;
        # follow-up issue to reconcile (Phase 2, #292).
    end
end

@testset "XYh1D — LocalMagnetization(:z) (OBC)" begin
    let m = XYh1D(; Jx=1.0, Jy=1.0, h=0.0), β = 1.0, N = 12
        mz = QAtlas.fetch(m, LocalMagnetization(:z), OBC(N); beta=β)
        @test length(mz) == N
        @test all(isfinite, mz)
    end
end

@testset "XYh1D — site-local observables (X, Y, EnergyLocal)" begin
    let m = XYh1D(; Jx=1.0, Jy=0.5, h=0.6), β = 1.0, N = 12
        mx = QAtlas.fetch(m, LocalMagnetization{:x}(), OBC(N); beta=β)
        my = QAtlas.fetch(m, LocalMagnetization(:y), OBC(N); beta=β)
        @test length(mx) == N && all(iszero, mx)
        @test length(my) == N && all(iszero, my)

        e = QAtlas.fetch(m, EnergyLocal(), OBC(N); beta=β)
        @test length(e) == N
        @test all(isfinite, e)
        # Sum of local energies should equal total Energy{:total}
        E_total = QAtlas.fetch(m, Energy{:total}(), OBC(N); beta=β)
        @test isapprox(sum(e), E_total; atol=1e-6)
    end
end
