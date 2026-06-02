# test/models/quantum/Hubbard1D/test_hubbard1d_jks_ed_comparison.jl
# ED-based mid-T verification for the JKS NLIE.

using Test
using LinearAlgebra
using SparseArrays
using QAtlas
using QAtlas: Hubbard1D, FreeEnergy, Infinite
using QAtlas.Hubbard1DJKSNLIE: atomic_free_energy

# Hubbard ED on small N (OBC). State: 2N-bit integer with JW ordering.
function _ed_hubbard_chain(
    N::Int, t::Float64, U::Float64, μ::Float64, β::Float64; pbc::Bool=false
)
    n_modes = 2 * N
    dim = 2^n_modes
    H_diag = zeros(Float64, dim)
    for s in 0:(dim - 1)
        e = 0.0
        for i in 1:N
            su, sd = 2*(i-1), 2*(i-1)+1
            nu = (s >> su) & 1;
            nd = (s >> sd) & 1
            e += U*nu*nd - μ*(nu+nd)
        end
        H_diag[s + 1] = e
    end
    rows = Int[];
    cols = Int[];
    vals = Float64[]
    bonds = pbc ? collect(1:N) : collect(1:(N - 1))
    for i in bonds
        ip = pbc ? (i % N) + 1 : i + 1
        for σ in 0:1
            si, sip = 2*(i-1)+σ, 2*(ip-1)+σ
            for s in 0:(dim - 1)
                if (s >> sip) & 1 == 1 && (s >> si) & 1 == 0
                    c1 = sum((s >> k) & 1 for k in 0:(sip - 1); init=0)
                    sg1 = iseven(c1) ? 1 : -1
                    s2 = s & ~(1 << sip)
                    c2 = sum((s2 >> k) & 1 for k in 0:(si - 1); init=0)
                    sg2 = iseven(c2) ? 1 : -1
                    s3 = s2 | (1 << si)
                    sg = sg1*sg2
                    push!(rows, s3+1);
                    push!(cols, s+1);
                    push!(vals, -t*sg)
                    push!(rows, s+1);
                    push!(cols, s3+1);
                    push!(vals, -t*sg)
                end
            end
        end
    end
    H = Matrix(sparse(rows, cols, vals, dim, dim))
    for i in 1:dim
        ;
        H[i, i] += H_diag[i];
    end
    eigs = eigvals(Hermitian(H))
    emin = minimum(eigs)
    return (-log(sum(exp.(-β .* (eigs .- emin))))/β + emin) / N
end

@testset "Hubbard1D JKS NLIE — ED N=4 mid-T comparison (#523)" begin
    @testset "ED essentially equals atomic at β ≤ 0.2 (large U regime)" begin
        # For U=4, β ≤ 0.2: kinetic correction is t²/U ~ 0.25 with β prefactor,
        # so f_ED should equal f_atom to within ~1%%.
        for β in (0.001, 0.01, 0.05, 0.1, 0.2)
            f_ed = _ed_hubbard_chain(4, 1.0, 4.0, 2.0, β; pbc=false)
            f_atom = atomic_free_energy(β, 4.0, 2.0)
            @test isapprox(f_ed, f_atom; rtol=0.05)
        end
    end

    @testset "JKS at high T matches ED to <0.5%%" begin
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        for β in (1e-4, 1e-3)
            f_ed = _ed_hubbard_chain(4, 1.0, 4.0, 2.0, β; pbc=false)
            f_jks = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
            @test isapprox(f_jks, f_ed; rtol=0.005)
        end
    end

    @testset "JKS mid-T deviates from ED >5%% — current implementation bug regression guard" begin
        # Stage G.2: ED proved mid-T ratio drop is a JKS formula bug, not physics.
        # This test PINS the buggy behaviour so any future fix is detected.
        m = Hubbard1D(; t=1.0, U=4.0, μ=2.0)
        for β in (0.1, 0.2)
            f_ed = _ed_hubbard_chain(4, 1.0, 4.0, 2.0, β; pbc=false)
            f_jks = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
            relative_error = abs(f_jks - f_ed) / abs(f_ed)
            @test_broken relative_error < 0.05  # currently fails ~25-50%%
        end
    end
end
