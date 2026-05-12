# ─────────────────────────────────────────────────────────────────────────────
# Meta-verification: generic Bloch builder vs hardcoded TB formulas
#
# For every tight-binding model in QAtlas/src/models/quantum/tightbinding/,
# verify that the hand-derived closed-form Bloch spectrum matches the
# generic `bloch_tb_spectrum` builder (which reads the unit cell
# definition from `Lattice2D.get_unit_cell` and constructs H(k) with
# no model-specific knowledge).
#
# If both agree, the hardcoded formulas are correct AND the generic
# builder is correct — mutual cross-validation.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

@testset "Generic Bloch TB vs hardcoded formulas" begin
    test_sizes = [(3, 3), (3, 4), (4, 4)]

    @testset "Honeycomb (Graphene)" begin
        for (Lx, Ly) in test_sizes
            λ_generic = bloch_tb_spectrum(Honeycomb, Lx, Ly, 1.0)
            λ_atlas = QAtlas.fetch(Graphene(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=1.0)
            @test λ_generic ≈ λ_atlas atol = 1e-10
        end
    end

    @testset "Kagome" begin
        for (Lx, Ly) in test_sizes
            λ_generic = bloch_tb_spectrum(Kagome, Lx, Ly, 1.0)
            λ_atlas = QAtlas.fetch(
                QAtlas.Kagome(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=1.0
            )
            @test λ_generic ≈ λ_atlas atol = 1e-10
        end
    end

    @testset "Lieb" begin
        for (Lx, Ly) in test_sizes
            λ_generic = bloch_tb_spectrum(Lieb, Lx, Ly, 1.0)
            λ_atlas = QAtlas.fetch(
                QAtlas.Lieb(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=1.0
            )
            @test λ_generic ≈ λ_atlas atol = 1e-10
        end
    end

    @testset "Triangular" begin
        for (Lx, Ly) in test_sizes
            λ_generic = bloch_tb_spectrum(Triangular, Lx, Ly, 1.0)
            λ_atlas = QAtlas.fetch(
                QAtlas.Triangular(), TightBindingSpectrum(); Lx=Lx, Ly=Ly, t=1.0
            )
            @test λ_generic ≈ λ_atlas atol = 1e-10
        end
    end

    @testset "Generic builder = real-space ED (honeycomb sanity)" begin
        # Triple cross-check: generic Bloch = hardcoded Bloch = real-space ED
        lat = build_lattice(Honeycomb, 3, 3)
        H = build_tight_binding(lat, 1.0)
        λ_real = sort(eigvals(Symmetric(H)))
        λ_generic = bloch_tb_spectrum(Honeycomb, 3, 3, 1.0)
        @test λ_real ≈ λ_generic atol = 1e-10
    end

    @testset "t scaling (generic builder)" begin
        for t in (0.5, 2.0, 3.7)
            λ1 = bloch_tb_spectrum(Kagome, 3, 3, 1.0)
            λt = bloch_tb_spectrum(Kagome, 3, 3, t)
            @test λt ≈ t .* λ1 rtol = 1e-12
        end
    end
end
