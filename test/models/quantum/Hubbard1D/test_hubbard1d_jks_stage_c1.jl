# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c1.jl
#
# Stage C.1 regression for the JKS auxiliary-function container (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSAuxFunctions,
    JKSContourGrid,
    init_atomic_limit!,
    init_atomic_limit,
    atomic_free_energy

@testset "Hubbard1D — JKS Stage C.1 aux container (#523)" begin
    @testset "Construction + zero init" begin
        aux = JKSAuxFunctions(64)
        @test length(aux) == 64
        @test length(aux.b) == 64
        @test length(aux.b_bar) == 64
        @test length(aux.c) == 64
        @test length(aux.c_bar) == 64
        @test all(iszero, aux.b)
        @test all(iszero, aux.c_bar)
    end

    @testset "Construction validates N > 0" begin
        @test_throws DomainError JKSAuxFunctions(0)
        @test_throws DomainError JKSAuxFunctions(-1)
    end

    @testset "copy is deep" begin
        aux = JKSAuxFunctions(8)
        aux.b[1] = 1.0 + 0im
        aux.c_bar[5] = 2.0 - 3im
        aux2 = copy(aux)
        @test aux2.b[1] == aux.b[1]
        @test aux2.c_bar[5] == aux.c_bar[5]
        # Mutating the copy must not affect the original.
        aux2.b[1] = 99.0 + 0im
        @test aux.b[1] == 1.0 + 0im
        @test aux2.b[1] == 99.0 + 0im
    end

    @testset "init_atomic_limit! fills all four arrays with constants" begin
        aux = JKSAuxFunctions(16)
        init_atomic_limit!(aux, 1.0, 4.0, 2.0)  # half-filling beta=1, U=4
        # Each array should be constant across x.
        for arr in (aux.b, aux.b_bar, aux.c, aux.c_bar)
            @test all(isapprox(v, arr[1]; atol=1e-14) for v in arr)
        end
        # At h = 0, c_const and c_bar_const coincide (z_up == z_down).
        @test isapprox(aux.c[1], aux.c_bar[1]; atol=1e-14)
        # b and b_bar are equal at h = 0 (initial-guess symmetry).
        @test isapprox(aux.b[1], aux.b_bar[1]; atol=1e-14)
    end

    @testset "init_atomic_limit (grid wrapper) returns fresh aux" begin
        g = JKSContourGrid(32, 2pi/3)
        aux = init_atomic_limit(g, 1.0, 4.0, 2.0)
        @test length(aux) == g.N
        @test all(isapprox(v, aux.b[1]; atol=1e-14) for v in aux.b)
    end

    @testset "init_atomic_limit! validates beta > 0" begin
        aux = JKSAuxFunctions(8)
        @test_throws DomainError init_atomic_limit!(aux, 0.0, 4.0, 2.0)
        @test_throws DomainError init_atomic_limit!(aux, -1.0, 4.0, 2.0)
    end

    @testset "init values reduce to atomic free energy via Boltzmann weights" begin
        # Cross-check: the Boltzmann weights used in init_atomic_limit! reproduce
        # Z_site = 1 + 2 exp(beta mu) + exp(beta(2 mu - U)) and hence
        # f_atomic. This guards against silent mis-typed weights.
        beta, U, mu = 1.0, 4.0, 2.0
        aux = init_atomic_limit(JKSContourGrid(4, 2pi/3), beta, U, mu)
        # Reverse-engineer Z_site from c_const = z_up / Z_site.
        z_up = exp(beta * mu)
        z_down = z_up
        z_empty = 1.0
        z_double = exp(beta * (2 * mu - U))
        Z_site = z_empty + z_up + z_down + z_double
        # Stage C.22c paper-precise init: PH symmetry gives c = exp(-βU/2),
        # not z_up/Z_site. The (z_up+z_down)/(z_empty+z_double) is still correct for b.
        @test isapprox(real(aux.c[1]), exp(-beta * U / 2); atol=1e-12)
        @test isapprox(real(aux.b[1]), (z_up + z_down) / (z_empty + z_double); atol=1e-12)
        # And the (independent) atomic_free_energy from Stage A.
        f_atomic = atomic_free_energy(beta, U, mu)
        f_check = -log(Z_site) / beta
        @test isapprox(f_atomic, f_check; atol=1e-12)
    end

    @testset "Magnetic field breaks c / c_bar symmetry" begin
        aux = JKSAuxFunctions(4)
        init_atomic_limit!(aux, 1.0, 4.0, 2.0; h=0.5)
        # h > 0 favours up-spin so c (up) > c_bar (down).
        @test real(aux.c[1]) > real(aux.c_bar[1])
    end
end
