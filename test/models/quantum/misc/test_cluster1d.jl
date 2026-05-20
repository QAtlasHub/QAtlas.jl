# ─────────────────────────────────────────────────────────────────────────────
# Test: Cluster1D (1D Z₂×Z₂ SPT cluster Hamiltonian).
#
# Verifies the two closed-form Phase-1 observables on a stabiliser model.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Cluster1D — Energy{:per_site} = -J (Phase 1)" begin
    @test QAtlas.fetch(Cluster1D(), Energy{:per_site}(), Infinite()) == -1.0
    @test QAtlas.fetch(Cluster1D(; J=2.5), Energy{:per_site}(), Infinite()) == -2.5
    # Linear scaling with J
    for J in (0.5, 1.0, 3.0, 7.0)
        @test QAtlas.fetch(Cluster1D(; J=J), Energy{:per_site}(), Infinite()) == -J
    end
end

@testset "Cluster1D — MassGap = 2J (Phase 1)" begin
    @test QAtlas.fetch(Cluster1D(), MassGap(), Infinite()) == 2.0
    @test QAtlas.fetch(Cluster1D(; J=3.0), MassGap(), Infinite()) == 6.0
end

@testset "Cluster1D — rejects J ≤ 0 (Phase 1)" begin
    @test_throws DomainError Cluster1D(; J=0.0)
    @test_throws DomainError Cluster1D(; J=-1.5)
    # Also via kwarg override on fetch
    m = Cluster1D(; J=1.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); J=0.0)
    @test_throws DomainError QAtlas.fetch(m, MassGap(), Infinite(); J=-2.0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "Cluster1D — verification cards" begin
    for J in (0.5, 1.0, 3.0)
        verify(
            Cluster1D(; J=J),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=(-J),
            agree_within=1e-12,
            refs=["Cluster-state Hamiltonian: e0 = -J (exact stabiliser ground state)"],
        )
        verify(
            Cluster1D(; J=J),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=2J,
            agree_within=1e-12,
            refs=["Cluster model gap = 2J"],
        )
    end
end

# ── additional verification cards (#381 batch) ─────────────────────────────
@testset "Cluster1D — Energy + MassGap closed forms (#381 batch)" begin
    # Every stabilizer K_i = σ^z_{i-1} σ^x_i σ^z_{i+1} contributes -J in
    # the cluster-state GS ⇒ E₀/N = -J. Single K_i flip costs 2J ⇒ Δ = 2J.
    for J in (0.5, 1.0, 2.0, 3.7)
        verify(
            Cluster1D(; J=J),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=-J,
            agree_within=1e-14,
            refs=["Briegel-Raussendorf 2001: cluster state is +1 eigenstate of every K_i ⇒ E₀/N = -J"],
        )
        verify(
            Cluster1D(; J=J),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=2 * J,
            agree_within=1e-14,
            refs=["Cluster Hamiltonian: single stabilizer flip K_i: +1 → -1 costs 2J ⇒ Δ = 2J"],
        )
    end
end

