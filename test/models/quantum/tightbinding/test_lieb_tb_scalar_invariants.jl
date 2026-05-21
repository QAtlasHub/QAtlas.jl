# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/tightbinding/test_lieb_tb_scalar_invariants.jl
#
# Scalar-invariant cross-validation of the Lieb TB Bloch closed form:
# (a) tr(H²) = 2 t² · n_NN_bonds = 8 t² Lx Ly, (b) max(spectrum) = 2√2 t
# at the Γ point, (c) real-space ED via Lattice2D.
#
# Restores a piece of the WHY-plane coverage removed by PR #449
# (deleted test/verification/tightbinding/test_lieb_tight_binding.jl).
# Reuses TightBindingChecksum + TightBindingMaxEnergy from PR #453.
#
# n_NN_bonds for Lieb PBC: 1 corner site (coord 4) + 2 edge sites
# (coord 2 each) per unit cell ⇒ (4+2+2)/2 · Lx·Ly = 4·Lx·Ly.
#
# Hubs added (both new — INVENTORY had ZERO Lieb entries before this PR):
#   Lieb/TightBindingChecksum/Infinite   (literature_value: 8 t² Lx Ly)
#   Lieb/TightBindingChecksum/Infinite   (ed_finite_size:   Lattice2D)
#   Lieb/TightBindingMaxEnergy/Infinite  (literature_value: 2√2 t at Γ)
#
# Pure verify(); branches off cards/honeycomb-tb-scalar-invariants.
# Refs #381; restores coverage lost in #449.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra
using Lattice2D: build_lattice
using Lattice2D: Lieb as L2D_Lieb

include(joinpath(@__DIR__, "..", "..", "..", "util", "tight_binding.jl"))

@testset "Lieb TB — TightBindingChecksum literature pin (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            verify(
                QAtlas.Lieb(),
                TightBindingChecksum(),
                Infinite();
                route=:literature_value,
                independent=8.0 * t^2 * Lx * Ly,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "Lieb TB tr(H²) = 2 t² · n_NN_bonds; (coord 4 + 2 + 2)/2 = 4 bonds per cell ⇒ tr(H²) = 8 t² Lx Ly",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end

@testset "Lieb TB — TightBindingChecksum vs real-space ED (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            lat = build_lattice(L2D_Lieb, Lx, Ly)
            H = build_tight_binding(lat, t)
            @assert norm(H - H') < 1e-12 "build_tight_binding returned non-Hermitian H at Lx=$(Lx), Ly=$(Ly), t=$(t)"
            ed_checksum = sum(abs2, eigvals(Symmetric(H)))
            verify(
                QAtlas.Lieb(),
                TightBindingChecksum(),
                Infinite();
                route=:ed_finite_size,
                independent=ed_checksum,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "ED black-box: Lattice2D.build_lattice(Lieb,Lx,Ly) + test/util/tight_binding.jl build_tight_binding → eigvals → Σ λ²",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end

@testset "Lieb TB — TightBindingMaxEnergy Γ-point literature pin (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            verify(
                QAtlas.Lieb(),
                TightBindingMaxEnergy(),
                Infinite();
                route=:literature_value,
                independent=2 * sqrt(2) * t,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "Lieb TB at Γ: E(k=0) = 2 |t| √(1 + 1) = 2√2 |t| (Lieb 1989; Tasaki 1998 — flat band at 0 plus ±E(k))",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end
