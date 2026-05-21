# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/tightbinding/test_triangular_tb_scalar_invariants.jl
#
# Scalar-invariant cross-validation of the Triangular TB Bloch closed form:
# (a) tr(H²) = 2 t² · n_NN_bonds = 6 t² Lx Ly, (b) Bloch checksum and
# max(spectrum) vs real-space ED via Lattice2D.
#
# Restores a piece of the WHY-plane coverage removed by PR #449
# (deleted test/verification/tightbinding/test_triangular_tight_binding.jl).
# Reuses TightBindingChecksum + TightBindingMaxEnergy from PR #453.
#
# Triangular is NOT bipartite. n_NN_bonds = 6 NN/site × Lx·Ly / 2 = 3·Lx·Ly.
# max(spectrum) reaches +3·t only when both Lx, Ly are multiples of 3
# (K-point on grid); otherwise it is strictly less, so no Lx,Ly-independent
# literature pin is possible — the MaxEnergy hub is exercised only via
# the ED route here.
#
# Hubs added (both new — INVENTORY had ZERO Triangular entries before this PR):
#   Triangular/TightBindingChecksum/Infinite   (literature_value: 6 t² Lx Ly)
#   Triangular/TightBindingChecksum/Infinite   (ed_finite_size:   Lattice2D)
#   Triangular/TightBindingMaxEnergy/Infinite  (ed_finite_size:   Lattice2D)
#
# Pure verify(); branches off cards/honeycomb-tb-scalar-invariants.
# Refs #381; restores coverage lost in #449.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra
using Lattice2D: build_lattice
using Lattice2D: Triangular as L2D_Triangular

include(joinpath(@__DIR__, "..", "..", "..", "util", "tight_binding.jl"))

@testset "Triangular TB — TightBindingChecksum literature pin (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            verify(
                QAtlas.Triangular(),
                TightBindingChecksum(),
                Infinite();
                route=:literature_value,
                independent=6.0 * t^2 * Lx * Ly,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "Triangular TB tr(H²) = 2 t² · n_NN_bonds; coord 6 per site ⇒ 3·Lx·Ly bonds, tr(H²) = 6 t² Lx Ly",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end

@testset "Triangular TB — TightBindingChecksum vs real-space ED (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            lat = build_lattice(L2D_Triangular, Lx, Ly)
            H = build_tight_binding(lat, t)
            @assert norm(H - H') < 1e-12 "build_tight_binding returned non-Hermitian H at Lx=$(Lx), Ly=$(Ly), t=$(t)"
            ed_checksum = sum(abs2, eigvals(Symmetric(H)))
            verify(
                QAtlas.Triangular(),
                TightBindingChecksum(),
                Infinite();
                route=:ed_finite_size,
                independent=ed_checksum,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "ED black-box: Lattice2D.build_lattice(Triangular,Lx,Ly) + test/util/tight_binding.jl build_tight_binding → eigvals → Σ λ²",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end

@testset "Triangular TB — TightBindingMaxEnergy vs real-space ED (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            lat = build_lattice(L2D_Triangular, Lx, Ly)
            H = build_tight_binding(lat, t)
            @assert norm(H - H') < 1e-12 "build_tight_binding returned non-Hermitian H at Lx=$(Lx), Ly=$(Ly), t=$(t)"
            ed_max = maximum(eigvals(Symmetric(H)))
            verify(
                QAtlas.Triangular(),
                TightBindingMaxEnergy(),
                Infinite();
                route=:ed_finite_size,
                independent=ed_max,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "ED black-box: Lattice2D real-space H, max(eigvals); Triangular K-point at (2π/3,2π/3) is on grid only when both Lx, Ly are multiples of 3",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end
