# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: XCube — subextensive GSD = 2^(2(Lx+Ly+Lz)-3).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "XCube — GSD = 2^(2(Lx+Ly+Lz)-3)" begin
    cases = (
        (2, 2, 2, big(2)^9),
        (3, 3, 3, big(2)^15),
        (2, 3, 4, big(2)^15),
        (5, 5, 5, big(2)^27),
    )
    for (Lx, Ly, Lz, expected) in cases
        @test QAtlas.fetch(XCube(), GroundStateDegeneracy(), PBC(); Lx=Lx, Ly=Ly, Lz=Lz) ==
            expected
    end
end

@testset "XCube — log_2 GSD scales linearly in (Lx+Ly+Lz)" begin
    L1 = QAtlas.fetch(XCube(), GroundStateDegeneracy(), PBC(); Lx=2, Ly=2, Lz=2)
    L2 = QAtlas.fetch(XCube(), GroundStateDegeneracy(), PBC(); Lx=3, Ly=3, Lz=3)
    # log_2 ratio = (2*9 - 3) - (2*6 - 3) = 15 - 9 = 6
    @test L2 == L1 * 2^6
end

@testset "XCube — DomainError on L_α < 2" begin
    @test_throws DomainError QAtlas.fetch(
        XCube(), GroundStateDegeneracy(), PBC(); Lx=1, Ly=2, Lz=2
    )
    @test_throws DomainError QAtlas.fetch(
        XCube(), GroundStateDegeneracy(), PBC(); Lx=2, Ly=0, Lz=2
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "XCube — verification cards" begin
    # X-cube fracton GSD on the 3-torus: log2(GSD) = 2(Lx+Ly+Lz) - 3.
    for (Lx, Ly, Lz) in ((2, 2, 2), (3, 3, 3))
        verify(
            XCube(),
            GroundStateDegeneracy(),
            PBC();
            route=:second_closed_form,
            fetch_kw=(; Lx=Lx, Ly=Ly, Lz=Lz),
            independent=2.0^(2 * (Lx + Ly + Lz) - 3),
            agree_within=0.5,
            refs=["X-cube fracton order: log2 GSD = 2(Lx+Ly+Lz) - 3"],
        )
    end
end
