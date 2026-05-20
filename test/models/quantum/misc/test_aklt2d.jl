using QAtlas, Test

@testset "AKLT2D — frustration-free Energy = 0 (Phase 1)" begin
    @test QAtlas.fetch(AKLT2D(), Energy{:per_site}(), Infinite()) == 0.0
    # J-independence (Hamiltonian is sum of non-negative projectors with
    # GS annihilation, so scaling J doesn't shift the ground state energy
    # from zero)
    for J in (0.5, 1.0, 2.5, 7.0)
        @test QAtlas.fetch(AKLT2D(; J=J), Energy{:per_site}(), Infinite()) == 0.0
    end
end

@testset "AKLT2D — rejects J ≤ 0 (Phase 1)" begin
    @test_throws DomainError AKLT2D(; J=0.0)
    @test_throws DomainError AKLT2D(; J=-1.0)
    m = AKLT2D(; J=1.0)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); J=0.0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "AKLT2D — verification cards" begin
    for J in (0.5, 1.0, 2.5)
        verify(
            AKLT2D(; J=J),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-12,
            refs=["Frustration-free: sum of non-negative projectors => e0 = 0"],
        )
    end
end
# ── additional verification cards (#381 batch 4) ─────────────────────────
@testset "AKLT2D — Energy frustration-free (#381 batch 4)" begin
    # 2D AKLT on honeycomb: H is a sum of non-negative spin-3 projectors,
    # VBS state is in the kernel of every term ⇒ frustration-free
    # GS with exact e₀ = 0 (AKLT 1988).
    for J in (0.5, 1.0, 2.0)
        verify(
            AKLT2D(; J=J),
            Energy(:per_site),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=0,
            refs=["AKLT 1988: 2D AKLT honeycomb is frustration-free ⇒ e₀ = 0 exactly"],
        )
    end
end
