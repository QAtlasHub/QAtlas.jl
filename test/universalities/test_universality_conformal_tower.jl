# test/universalities/test_universality_conformal_tower.jl

using QAtlas, Test

@testset "Universality: CFT Conformal Tower of States" begin
    # ─── Ising PBC ────────────────────────────────────────────────────────────
    @testset "Ising PBC" begin
        verify(
            Universality(:Ising),
            ConformalTower(),
            PBC();
            route=:second_closed_form,
            independent=(0.25 * π * 2.0 / 16.0),
            agree_within=1e-12,
            refs=[
                "Cardy 1986: Ising PBC spin state dimension Δ = 1/8",
            ],
            fetch_kw=(; L=16.0, v=2.0),
            subject_extract=raw -> raw[2].energy,
        )

        tower = QAtlas.fetch(Universality(:Ising), ConformalTower(), PBC(); L=16.0, v=2.0)
        @test length(tower) == 4
        @test tower[1] == (energy=0.0, dimension=0.0, degeneracy=1)
        @test tower[2].dimension == 0.125
        @test tower[2].degeneracy == 1
        @test tower[3].dimension == 1.0
        @test tower[3].degeneracy == 1
        @test tower[4].dimension == 2.0
        @test tower[4].degeneracy == 2
    end

    # ─── Ising OBC ────────────────────────────────────────────────────────────
    @testset "Ising OBC" begin
        verify(
            Universality(:Ising),
            ConformalTower(),
            OBC();
            route=:second_closed_form,
            independent=(π * 2.0 / 256.0), # = 0.0625 * π * v / L = π * 2 / (16 * 16)
            agree_within=1e-12,
            refs=[
                "Cardy 1986 / Blöte-Cardy-Nightingale 1986: Ising OBC free-free spin state dimension h = 1/16",
            ],
            fetch_kw=(; L=16.0, v=2.0),
            subject_extract=raw -> raw[2].energy,
        )

        tower = QAtlas.fetch(Universality(:Ising), ConformalTower(), OBC(); L=16.0, v=2.0)
        @test length(tower) == 3
        @test tower[1] == (energy=0.0, dimension=0.0, degeneracy=1)
        @test tower[2].dimension == 0.0625
        @test tower[2].degeneracy == 1
        @test tower[3].dimension == 0.5
        @test tower[3].degeneracy == 1
    end

    # ─── Heisenberg PBC ───────────────────────────────────────────────────────
    @testset "Heisenberg PBC" begin
        verify(
            Universality(:Heisenberg),
            ConformalTower(),
            PBC();
            route=:second_closed_form,
            independent=(0.5 * π * 2.0 / 16.0), # = 0.25 * 2π v / L
            agree_within=1e-12,
            refs=[
                "Cardy 1986 / Affleck 1986: Heisenberg PBC spinon state dimension Δ = 1/4",
            ],
            fetch_kw=(; L=16.0, v=2.0),
            subject_extract=raw -> raw[2].energy,
        )

        tower = QAtlas.fetch(
            Universality(:Heisenberg), ConformalTower(), PBC(); L=16.0, v=2.0
        )
        @test length(tower) == 3
        @test tower[1] == (energy=0.0, dimension=0.0, degeneracy=1)
        @test tower[2].dimension == 0.25
        @test tower[2].degeneracy == 4
        @test tower[3].dimension == 1.0
        @test tower[3].degeneracy == 9
    end

    # ─── Error path and domain validation ─────────────────────────────────────
    @testset "Error and validation paths" begin
        # Invalid arguments L, v
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), ConformalTower(), PBC(); L=-5.0, v=2.0
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), ConformalTower(), PBC(); L=16.0, v=0.0
        )

        # Unsupported classes
        @test_throws ErrorException QAtlas.fetch(
            Universality(:XY), ConformalTower(), PBC(); L=16.0, v=2.0
        )

        # Heisenberg OBC not supported
        @test_throws ErrorException QAtlas.fetch(
            Universality(:Heisenberg), ConformalTower(), OBC(); L=16.0, v=2.0
        )
    end
end

@testset "Model: Conformal Tower of States Delegates" begin
    # ─── TFIM ─────────────────────────────────────────────────────────────────
    @testset "TFIM" begin
        tfim_crit = TFIM(; J=1.5, h=1.5)
        # PBC
        tower_pbc = QAtlas.fetch(tfim_crit, ConformalTower(), PBC(16))
        expected_v = 3.0 # 2J
        scale_pbc = 2 * π * expected_v / 16
        @test tower_pbc[2].energy ≈ 0.125 * scale_pbc
        @test tower_pbc[2].dimension == 0.125

        # OBC
        tower_obc = QAtlas.fetch(tfim_crit, ConformalTower(), OBC(16))
        scale_obc = π * expected_v / 16
        @test tower_obc[2].energy ≈ 0.0625 * scale_obc
        @test tower_obc[2].dimension == 0.0625

        # Off-critical throws DomainError
        tfim_off = TFIM(; J=1.5, h=1.0)
        @test_throws DomainError QAtlas.fetch(tfim_off, ConformalTower(), PBC(16))
    end

    # ─── Heisenberg1D ─────────────────────────────────────────────────────────
    @testset "Heisenberg1D" begin
        heis = Heisenberg1D()
        tower = QAtlas.fetch(heis, ConformalTower(), PBC(16); J=1.5)
        expected_v = 1.5 * π / 2
        scale = 2 * π * expected_v / 16
        @test tower[2].energy ≈ 0.25 * scale
        @test tower[2].dimension == 0.25

        # OBC throws ErrorException
        @test_throws ErrorException QAtlas.fetch(heis, ConformalTower(), OBC(16); J=1.5)
    end
end
