# ─────────────────────────────────────────────────────────────────────────────
# AKLT1D — β = ∞ (T = 0) thermodynamic-limit fetches
#
# The AKLT model is not Bethe-ansatz integrable, so no closed-form
# analytic / TM / integral-equation reduction is known for the finite-T
# thermodynamic functions (free energy, entropy, specific heat,
# susceptibility) at arbitrary β.  Only the β = ∞ ground-state limits
# are registered in src/, all derived from the AKLT bond-projector
# decomposition (AKLT 1988, doi:10.1007/BF01218021).
#
# This file:
#   1. checks the registered β = ∞ values against the analytic formulas,
#   2. cross-checks OBC against an independent dense-ED of the OBC
#      Hamiltonian (the 4-fold edge-mode degeneracy is also verified),
#   3. confirms finite-β fetches throw DomainError everywhere,
#   4. confirms OBC × SusceptibilityZZ throws at every β (Curie tail
#      divergence at β = ∞).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra
using LinearAlgebra: Hermitian, eigvals

@testset "AKLT1D β=∞ — OBC FreeEnergy = -(2J/3)(N-1)/N" begin
    # ── Analytic identity vs registered fetch ───────────────────────────────
    for J in (0.5, 1.0, 2.5), N in (2, 3, 4, 5, 6, 8)
        f_fetch = QAtlas.fetch(AKLT1D(; J=J), FreeEnergy(), OBC(N); beta=Inf)
        f_closed = -(2.0 / 3.0) * J * (N - 1) / N
        @test f_fetch ≈ f_closed atol = 1e-14
    end

    # ── Independent dense-ED cross-check (bond-projector says E_GS_OBC =
    #    -(2J/3)(N-1) exactly; ED on the OBC Hamiltonian should agree). ──────
    for J in (0.5, 1.0), N in (3, 4, 5, 6)
        H = QAtlas._aklt_hamiltonian_matrix(AKLT1D(; J=J), N, OBC(N))
        E_gs = real(eigvals(Hermitian(H))[1])
        f_ed = E_gs / N
        f_fetch = QAtlas.fetch(AKLT1D(; J=J), FreeEnergy(), OBC(N); beta=Inf)
        @test f_fetch ≈ f_ed atol = 1e-10
    end
end

@testset "AKLT1D β=∞ — OBC ThermalEntropy = log(4)/N (edge degeneracy)" begin
    for N in (2, 3, 4, 5, 6, 8)
        s = QAtlas.fetch(AKLT1D(), ThermalEntropy(), OBC(N); beta=Inf)
        @test s ≈ log(4) / N atol = 1e-14
    end

    # ── Confirm the 4-fold edge-mode GS degeneracy that underlies log(4) ───
    # The next state above the 4-fold manifold should be separated by a
    # finite N-dependent gap (the OBC AKLT spectrum: 4 (nearly) degenerate
    # ground states from spin-½ edge modes, then a Haldane-like gap).
    for N in (3, 4, 5, 6)
        H = QAtlas._aklt_hamiltonian_matrix(AKLT1D(), N, OBC(N))
        ev = sort(real.(eigvals(Hermitian(H))))
        @test ev[4] - ev[1] < 1e-9                # 4-fold manifold
        @test ev[5] - ev[1] > 0.1                 # gapped to next state
    end
end

@testset "AKLT1D β=∞ — OBC SpecificHeat = 0" begin
    for J in (0.5, 1.0, 2.5), N in (2, 3, 5, 6)
        c = QAtlas.fetch(AKLT1D(; J=J), SpecificHeat(), OBC(N); beta=Inf)
        @test c == 0.0
    end
end

@testset "AKLT1D β=∞ — PBC FreeEnergy = -2J/3 (N-independent)" begin
    for J in (0.5, 1.0, 2.5), N in (2, 3, 4, 5, 6, 8)
        f = QAtlas.fetch(AKLT1D(; J=J), FreeEnergy(), PBC(N); beta=Inf)
        @test f ≈ -(2.0 / 3.0) * J atol = 1e-14
    end
end

@testset "AKLT1D β=∞ — PBC entropy / heat / susceptibility = 0" begin
    for N in (2, 3, 4, 5)
        @test QAtlas.fetch(AKLT1D(), ThermalEntropy(), PBC(N); beta=Inf) == 0.0
        @test QAtlas.fetch(AKLT1D(), SpecificHeat(), PBC(N); beta=Inf) == 0.0
        @test QAtlas.fetch(AKLT1D(), SusceptibilityZZ(), PBC(N); beta=Inf) == 0.0
    end
end

@testset "AKLT1D β=∞ — Infinite limits match Infinite GSED" begin
    for J in (0.5, 0.7, 1.0, 1.3, 2.5)
        m = AKLT1D(; J=J)
        f_inf = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=Inf)
        e_gs = QAtlas.fetch(m, GroundStateEnergyDensity(), Infinite())
        @test f_inf ≈ e_gs atol = 1e-14
        @test f_inf ≈ -(2.0 / 3.0) * J atol = 1e-14
        @test QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=Inf) == 0.0
        @test QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=Inf) == 0.0
        @test QAtlas.fetch(m, SusceptibilityZZ(), Infinite(); beta=Inf) == 0.0
    end
end

@testset "AKLT1D — finite β throws DomainError (no analytic reduction)" begin
    m = AKLT1D()
    for bc in (OBC(4), PBC(4), Infinite())
        for q in (FreeEnergy(), ThermalEntropy(), SpecificHeat())
            for beta in (1.0, 0.5, 2.0, 10.0, 0.0)
                @test_throws DomainError QAtlas.fetch(m, q, bc; beta=beta)
            end
        end
    end
    # SusceptibilityZZ: only PBC and Infinite return 0 at β=Inf; finite β throws.
    for bc in (PBC(4), Infinite())
        for beta in (1.0, 0.5, 2.0, 0.0)
            @test_throws DomainError QAtlas.fetch(m, SusceptibilityZZ(), bc; beta=beta)
        end
    end
end

@testset "AKLT1D — OBC SusceptibilityZZ always throws (Curie divergence)" begin
    # The 4-fold OBC ground manifold gives ⟨(S^z_tot)²⟩_GS = 1/2 (one
    # singlet at 0 + three triplet states at 1, 0, 1, averaged), so
    # χ_OBC(β → ∞) ≈ β / (2N) diverges.  fetch refuses to return ∞.
    m = AKLT1D()
    for beta in (Inf, 100.0, 10.0, 1.0, 0.5, 0.0), N in (2, 3, 4)
        @test_throws DomainError QAtlas.fetch(m, SusceptibilityZZ(), OBC(N); beta=beta)
    end
end

# ── verify() cards (registry round-trip) ───────────────────────────────────
# These exercise the route=:second_closed_form pathway against the
# documented analytic forms, so a registry mismatch (wrong reliability,
# broken @register row, hub-dispatch regression) fails the suite alongside
# the @testset block above.
@testset "AKLT1D β=∞ — verification cards" begin
    for J in (0.7, 1.0, 1.3)
        verify(
            AKLT1D(; J=J),
            FreeEnergy(),
            OBC(4);
            route=:second_closed_form,
            fetch_kw=(; beta=Inf),
            independent=-(2.0 / 3.0) * J * 3 / 4,    # N=4: -(2J/3)·3/4
            agree_within=1e-14,
            refs=["AKLT 1988: bond-projector E_GS_OBC = -(2J/3)(N-1)"],
        )
        verify(
            AKLT1D(; J=J),
            FreeEnergy(),
            PBC(4);
            route=:second_closed_form,
            fetch_kw=(; beta=Inf),
            independent=-(2.0 / 3.0) * J,
            agree_within=1e-14,
            refs=["AKLT 1988: PBC VBS unique, f = -2J/3"],
        )
        verify(
            AKLT1D(; J=J),
            FreeEnergy(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=Inf),
            independent=-(2.0 / 3.0) * J,
            agree_within=1e-14,
            refs=["AKLT 1988: bulk f(∞) = -2J/3"],
        )
    end
    for N in (3, 4, 6)
        verify(
            AKLT1D(),
            ThermalEntropy(),
            OBC(N);
            route=:second_closed_form,
            fetch_kw=(; beta=Inf),
            independent=log(4) / N,
            agree_within=1e-14,
            refs=["AKLT 1988: 4-fold edge-mode GS degeneracy, s(∞) = log 4 / N"],
        )
    end
end
