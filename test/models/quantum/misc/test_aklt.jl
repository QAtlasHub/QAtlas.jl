# ─────────────────────────────────────────────────────────────────────────────
# Test: AKLT1D — exact VBS ground state + Haldane gap
#
# Values are verified by the verify() cards below (the new system). This
# file retains ONLY the structural / error / identity / relational guards
# that verify() architecturally cannot express (constructor invariants,
# OBC spectrum shape & 4-fold edge-state degeneracy, ZZ-correlation
# symmetry/sign/ratio identities, structure-factor periodicity & FT-of-
# correlation identity). Legacy hand-rolled value @testsets are deleted
# — superseded by the verification cards.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "AKLT1D — structural / OBC-spectrum / identity guards" begin
    @testset "Constructor variants" begin
        @test AKLT1D().J == 1.0
        @test AKLT1D(; J=2.5).J == 2.5
    end

    @testset "OBC ExactSpectrum shape (sorted, real, length 3^N)" begin
        m = AKLT1D(; J=1.0)
        for N in (2, 3, 4)
            λ = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
            @test length(λ) == 3^N
            @test issorted(λ)
            @test all(isreal, λ)
        end
    end

    @testset "4-fold OBC ground-state degeneracy (S_tot = 0 ⊕ 1)" begin
        # AKLT theorem: two free spin-1/2 edge modes → 4 = 1 + 3
        # degenerate ground states. Structural; not a single value.
        m = AKLT1D(; J=1.0)
        for N in (4, 6, 8)
            λ = QAtlas.fetch(m, ExactSpectrum(), OBC(N))
            @test λ[4] - λ[1] < 1e-8                   # 4-fold manifold
            @test λ[5] - λ[4] > 0.05                    # well separated
        end
    end

    @testset "OBC ground-state energy scaling (analytic + bounds)" begin
        # E_0(N) = -(2/3)(N-1) J exact for OBC AKLT (every bond projector
        # gives zero). e_0/N → -2/3 with a 1/N edge correction.
        m = AKLT1D(; J=1.0)
        e_inf = -2 / 3
        for N in (4, 6, 8)
            E0 = QAtlas.fetch(m, ExactSpectrum(), OBC(N))[1]
            @test E0 ≈ -(2 / 3) * (N - 1) atol = 1e-10
            @test E0 / N > e_inf
            @test abs(E0 / N - e_inf) ≤ 2 / (3N)
        end
    end

    @testset "ZZCorrelation — symmetry, sign-alternation, ratio, J-independence" begin
        # Relational identities (not single-value cards):
        # r ↔ -r symmetry; alternating sign; ratio c_{r+1}/c_r = -1/3;
        # J-independent (VBS ground state J-invariant for J>0).
        m = AKLT1D(; J=1.0)
        for r in 1:5
            @test QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=(-r)) ≈
                QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=r) atol = 1e-14
        end
        c1 = QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=1)
        c2 = QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=2)
        c3 = QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=3)
        @test c1 < 0 && c2 > 0 && c3 < 0
        @test abs(c2 / c1) ≈ 1 / 3 atol = 1e-14
        @test abs(c3 / c2) ≈ 1 / 3 atol = 1e-14
        for J in (0.3, 1.0, 4.2)
            @test QAtlas.fetch(
                AKLT1D(; J=J), ZZCorrelation(; mode=:static), Infinite(); r=2
            ) ≈ QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=2) atol = 1e-14
        end
    end

    @testset "ZZStructureFactor — monotonicity, periodicity, J-independence" begin
        # Relational: 0 < S(q) < S(π) for 0 < q < π; 2π-periodic; even;
        # J-independent.
        m = AKLT1D(; J=1.0)
        Sπ = QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=π)
        for q in range(0.1, π - 0.1; length=12)
            S = QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=q)
            @test 0.0 < S < 2.0
            @test S < Sπ
        end
        @test QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=0.7) ≈
            QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=0.7 + 2π) atol = 1e-12
        @test QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=-0.7) ≈
            QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=0.7) atol = 1e-14
        @test QAtlas.fetch(AKLT1D(; J=2.6), ZZStructureFactor(), Infinite(); q=π) ≈ 2.0 atol =
            1e-14
    end

    @testset "Structure factor IS the Fourier transform of the correlation" begin
        # Cross-quantity identity: S_zz(q) = Σ_r e^{iqr} ⟨Sᶻ₀Sᶻ_r⟩.
        # Truncated lattice sum (geometric tail < 1e-12 by r = 40).
        m = AKLT1D(; J=1.0)
        for q in (0.3, 1.0, 2.0, π)
            Ssum = QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=0)
            for r in 1:40
                Ssum +=
                    2 *
                    cos(q * r) *
                    QAtlas.fetch(m, ZZCorrelation(; mode=:static), Infinite(); r=r)
            end
            @test Ssum ≈ QAtlas.fetch(m, ZZStructureFactor(), Infinite(); q=q) atol = 1e-9
        end
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
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
@testset "AKLT1D — closed-form / DMRG cards (#381 batch)" begin
    # GroundStateEnergyDensity/Infinite: AKLT 1988 exact VBS GS energy
    # density e₀ = -2J/3, J-linear, J-independent of the wavefunction.
    for J in (0.5, 1.0, 2.0, 3.7)
        verify(
            AKLT1D(; J=J),
            GroundStateEnergyDensity(),
            Infinite();
            route=:second_closed_form,
            independent=-2 * J / 3,
            agree_within=1e-14,
            refs=["AKLT 1988: VBS ground state is the exact null space of every bond P₂ projector ⇒ e₀ = -2J/3"],
        )
    end

    # CorrelationLength/Infinite: ξ = 1/log 3 (AKLT 1988); J-independent
    # because the VBS wavefunction does not depend on J > 0.
    for J in (0.5, 1.0, 2.0)
        verify(
            AKLT1D(; J=J),
            CorrelationLength(),
            Infinite();
            route=:second_closed_form,
            independent=1 / log(3),
            agree_within=1e-14,
            refs=["AKLT 1988: ⟨S^z_0 S^z_r⟩ = (-1)^r (4/3) 3^{-r} ⇒ ξ = 1/log 3"],
        )
    end

    # StringOrderParameter/Infinite: O_str = 4/9 (AKLT 1988, Kennedy-Tasaki
    # 1992 hidden Z2×Z2 symmetry-breaking order parameter); J-independent.
    for J in (0.5, 1.0, 2.0)
        verify(
            AKLT1D(; J=J),
            StringOrderParameter(),
            Infinite();
            route=:second_closed_form,
            independent=4 / 9,
            agree_within=1e-14,
            refs=["Kennedy-Tasaki 1992 on AKLT VBS: O_str = -⟨S^z_i e^{iπ Σ S^z_k} S^z_j⟩ → 4/9 at r → ∞"],
        )
    end

    # MassGap/Infinite: Haldane gap Δ ≈ 0.350 J — DMRG literature value
    # (García-Saez/Murg/Verstraete 2013, PRB 88, 245118); no closed form.
    #
    # IMPORTANT — agree_within=5e-3 is the DMRG literature uncertainty floor
    # (García-Saez–Murg–Verstraete report Δ/J ≈ 0.3502 ± 0.0001; ~25× their
    # quoted uncertainty gives headroom for hub-stored precision choices
    # between 0.350 and 0.35048). This tolerance is INTENTIONALLY loose and
    # should NOT be tightened by future maintainers without a new high-
    # precision DMRG/MERA result superseding GMV 2013.
    for J in (1.0, 2.0)
        verify(
            AKLT1D(; J=J),
            MassGap(),
            Infinite();
            route=:literature_value,
            independent=0.350 * J,
            agree_within=5e-3,
            refs=["García-Saez–Murg–Verstraete 2013 (PRB 88, 245118): AKLT Haldane gap Δ ≈ 0.350 J (DMRG)"],
        )
    end
end

