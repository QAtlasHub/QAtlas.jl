# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/tightbinding/test_kagome_tb_scalar_invariants.jl
#
# Scalar-invariant cross-validation of the Kagome TB Bloch closed form
# against (a) tr(H²) = 2t²·n_NN_bonds = 12 t² Lx Ly and (b) the kagome
# flat-band value max(spectrum) = +2t, plus (c) real-space ED.
#
# Restores a piece of the WHY-plane coverage removed by PR #449:
# the deleted test/verification/tightbinding/test_kagome_tight_binding.jl
# cross-checked the Kagome closed-form Bloch spectrum against real-
# space ED. Sister to PR #453 (Honeycomb); reuses TightBindingChecksum
# and TightBindingMaxEnergy quantities defined there.
#
# n_NN_bonds for kagome PBC: each unit cell has 3 sites with coordination
# 4 (two within-cell + two cross-cell), giving 6·Lx·Ly unique NN bonds.
#
# Hubs added (both new — INVENTORY had ZERO Kagome entries before this PR):
#   Kagome/TightBindingChecksum/Infinite  (literature_value: 12 t² Lx Ly)
#   Kagome/TightBindingChecksum/Infinite  (ed_finite_size:   Lattice2D)
#   Kagome/TightBindingMaxEnergy/Infinite (literature_value: 2 t — flat band)
#
# Pure verify(); branches off cards/honeycomb-tb-scalar-invariants
# (which adds the TightBindingChecksum + TightBindingMaxEnergy quantities).
# Refs #381; restores coverage lost in #449.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra
using Lattice2D: build_lattice
using Lattice2D: Kagome as L2D_Kagome

include(joinpath(@__DIR__, "..", "..", "..", "util", "tight_binding.jl"))

@testset "Kagome TB — TightBindingChecksum literature pin (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            verify(
                QAtlas.Kagome(),
                TightBindingChecksum(),
                Infinite();
                route=:literature_value,
                independent=12.0 * t^2 * Lx * Ly,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "Kagome TB tr(H²) = 2 t² · n_NN_bonds; coord=4 per site × 3 sites/cell ÷ 2 ⇒ 6·Lx·Ly bonds, tr(H²) = 12 t² Lx Ly",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end

@testset "Kagome TB — TightBindingChecksum vs real-space ED (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            lat = build_lattice(L2D_Kagome, Lx, Ly)
            H = build_tight_binding(lat, t)
            # Hermitian guard: see PR #453 A3.
            @assert norm(H - H') < 1e-12 "build_tight_binding returned non-Hermitian H at Lx=$(Lx), Ly=$(Ly), t=$(t)"
            ed_checksum = sum(abs2, eigvals(Symmetric(H)))
            verify(
                QAtlas.Kagome(),
                TightBindingChecksum(),
                Infinite();
                route=:ed_finite_size,
                independent=ed_checksum,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "ED black-box: Lattice2D.build_lattice(Kagome,Lx,Ly) + test/util/tight_binding.jl build_tight_binding → eigvals → Σ λ²",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end

@testset "Kagome TB — TightBindingMaxEnergy flat-band literature pin (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            verify(
                QAtlas.Kagome(),
                TightBindingMaxEnergy(),
                Infinite();
                route=:literature_value,
                independent=2.0 * t,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "Kagome TB exactly flat band at +2t (Syôzi 1951; Bergman-Wu-Balents 2008): max(spectrum) = 2·|t|, independent of Lx, Ly",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end
