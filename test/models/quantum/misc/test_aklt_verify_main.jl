# ─────────────────────────────────────────────────────────────────────────────
# AKLT1D — verification cards (split, formerly testset 2 of 3)
#
# Split out of test/models/quantum/misc/test_aklt.jl (5.9 min on s02) so
# the three top-level testsets each run on their own shard. Helpers
# spin_ops, chain_hamiltonian, two_point, verify_profile_Ns come from
# test/util/{generic_ed,verify}.jl via runtests.jl ambient include.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

@testset "AKLT1D — verification cards" begin
    using LinearAlgebra: Hermitian, eigen, kron

    # Spin-1 site operators + nearest-neighbour AKLT bond term
    # H_b = S·S + (1/3)(S·S)^2.  Black-box ED via chain_hamiltonian, never
    # a QAtlas internal builder.
    Sx, Sy, Sz = spin_ops(1)
    SS = kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz)
    bond = SS + (1 / 3) * (SS * SS)
    Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8), nightly=(6, 8))

    function aklt_manifold(N)
        F = eigen(Hermitian(Matrix(chain_hamiltonian(3, N, bond))))
        E0 = F.values[1]
        deg = findall(e -> isapprox(e, E0; atol=1e-8), F.values)
        return F.values, [F.vectors[:, k] for k in deg]
    end
    function avg_zz(vecs, N, i, j)
        return sum(two_point(psi, 3, N, Sz, i, j) for psi in vecs) / length(vecs)
    end

    # ── closed-form analytical values (replace legacy J-loop value tests) ──
    # GSED Infinite: closed form -2J/3 (bond-projector decomposition).
    for J in (0.5, 1.0, 2.5)
        verify(
            AKLT1D(; J=J),
            GroundStateEnergyDensity(),
            Infinite();
            route=:second_closed_form,
            independent=-2 * J / 3,
            agree_within=1e-14,
            refs=["Affleck-Kennedy-Lieb-Tasaki 1988: e0 = -2J/3"],
        )
    end

    # CorrelationLength Infinite: ξ = 1/log 3 (VBS transfer matrix), J-indep.
    verify(
        AKLT1D(),
        CorrelationLength(),
        Infinite();
        route=:second_closed_form,
        independent=1 / log(3),
        agree_within=1e-12,
        refs=["AKLT 1988: ξ = 1/log 3 from VBS transfer matrix"],
    )

    # StringOrderParameter Infinite: O_str = 4/9 (Kennedy-Tasaki 1992), J-indep.
    verify(
        AKLT1D(),
        StringOrderParameter(),
        Infinite();
        route=:second_closed_form,
        independent=4 / 9,
        agree_within=1e-14,
        refs=["Kennedy-Tasaki 1992: O_str = 4/9 (hidden Z₂×Z₂ in Haldane phase)"],
    )

    # MassGap Infinite: García-Saez-Murg-Verstraete 2013 DMRG ≈ 0.350 J.
    for J in (0.5, 1.0, 2.5)
        verify(
            AKLT1D(; J=J),
            MassGap(),
            Infinite();
            route=:literature_value,
            independent=0.350 * J,
            agree_within=1e-3 * max(1.0, J),
            refs=["García-Saez-Murg-Verstraete 2013: Haldane gap Δ ≈ 0.350 J (DMRG)"],
        )
    end

    # ZZCorrelation Infinite: AKLT-1988 closed form (-1)^r (4/3) 3^{-|r|}
    # (and 2/3 at r=0). r=2 also cross-checked by ed_finite_size below.
    for r in (0, 1, 2, 3, 4)
        closed = r == 0 ? 2 / 3 : (-1)^r * (4 / 3) * 3.0^(-r)
        verify(
            AKLT1D(),
            ZZCorrelation(; mode=:static),
            Infinite();
            route=:second_closed_form,
            independent=closed,
            agree_within=1e-14,
            refs=["AKLT 1988: ⟨Sᶻ₀Sᶻ_r⟩ = (-1)^r (4/3) 3^{-|r|}, J-independent"],
            fetch_kw=(; r=r),
        )
    end

    # ZZStructureFactor Infinite: closed form S(q) = 2(1-cos q)/(5+3cos q).
    for q in (0.0, π / 2, π)
        closed = q == 0.0 ? 0.0 : 2 * (1 - cos(q)) / (5 + 3 * cos(q))
        verify(
            AKLT1D(),
            ZZStructureFactor(),
            Infinite();
            route=:second_closed_form,
            independent=closed,
            agree_within=1e-14,
            refs=["Arovas-Auerbach-Haldane 1988: S_zz(q) = 2(1-cos q)/(5+3cos q)"],
            fetch_kw=(; q=q),
        )
    end

    # ── independent finite-size dense-ED cards (non-circular vs closed forms) ──
    verify(
        AKLT1D(),
        Energy(:per_site),
        Infinite();
        route=:ed_finite_size,
        independent=[aklt_manifold(N)[1][1] / (N - 1) for N in Ns],
        at=["N=$N" for N in Ns],
        agree_within=1e-9,
        refs=["Affleck-Kennedy-Lieb-Tasaki 1988"],
    )

    let r = 2
        ind = Float64[]
        for N in Ns
            _, vecs = aklt_manifold(N)
            i0 = cld(N, 2)
            push!(ind, avg_zz(vecs, N, i0, i0 + r))
        end
        verify(
            AKLT1D(),
            ZZCorrelation(; mode=:static),
            Infinite();
            route=:ed_finite_size,
            independent=ind,
            at=["N=$N" for N in Ns],
            agree_within=3e-2,
            fetch_kw=(; r=r),
            refs=["Affleck-Kennedy-Lieb-Tasaki 1988"],
        )
    end

    let q = pi
        ind = Float64[]
        for N in Ns
            _, vecs = aklt_manifold(N)
            i0 = cld(N, 2)
            S = avg_zz(vecs, N, i0, i0)
            for r in 1:(N - i0)
                S += 2 * cos(q * r) * avg_zz(vecs, N, i0, i0 + r)
            end
            push!(ind, S)
        end
        verify(
            AKLT1D(),
            ZZStructureFactor(),
            Infinite();
            route=:ed_finite_size,
            independent=ind,
            at=["N=$N" for N in Ns],
            agree_within=0.15,
            fetch_kw=(; q=q),
            refs=["Arovas-Auerbach-Haldane 1988"],
        )
    end
end

# ── additional verification cards (#381 batch) ─────────────────────────────
