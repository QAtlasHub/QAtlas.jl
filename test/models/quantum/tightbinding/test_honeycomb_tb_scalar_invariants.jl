# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/tightbinding/test_honeycomb_tb_scalar_invariants.jl
#
# Scalar-invariant cross-validation of the Honeycomb TB Bloch closed form
# against (a) the analytical chiral-symmetry identity tr(H²) = 2t²·n_NN
# and (b) real-space ED via Lattice2D.build_tight_binding.
#
# Restores a piece of the WHY-plane coverage removed by PR #449 phase 1:
# the deleted test/verification/tightbinding/test_graphene_tight_binding.jl
# cross-checked the Honeycomb closed-form Bloch spectrum against real-
# space ED for several (Lx, Ly) sizes. The full Bloch spectrum is a
# Vector, so verify() cannot pin it directly — this file pins the scalar
# reduction TightBindingChecksum = Σ λᵢ² = tr(H²), which is (1) sensitive
# to typos in the Bloch closed form (any sign or factor slip changes the
# checksum), (2) exact at machine precision against the literature
# identity 2t²·n_NN_bonds for any chiral bipartite lattice, and (3) also
# exact against real-space ED to the same precision.
#
# n_NN_bonds for honeycomb PBC: each unit cell has 1 A site with 3 NN
# bonds to B sites, never double-counted ⇒ n_NN_bonds = 3·Lx·Ly.
#
# Hubs added (both new — INVENTORY had ZERO Honeycomb entries before this PR):
#   Honeycomb/TightBindingChecksum/Infinite (literature_value: 2t²·3·Lx·Ly)
#   Honeycomb/TightBindingChecksum/Infinite (ed_finite_size: Lattice2D)
#
# Pure verify(); branches off main. Refs #381; restores coverage lost in #449.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra
using Lattice2D: build_lattice
using Lattice2D: Honeycomb as L2D_Honeycomb

include(joinpath(@__DIR__, "..", "..", "..", "util", "tight_binding.jl"))

@testset "Honeycomb TB — TightBindingChecksum literature pin (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            verify(
                QAtlas.Honeycomb(),
                TightBindingChecksum(),
                Infinite();
                route=:literature_value,
                independent=2.0 * t^2 * 3 * Lx * Ly,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "Chiral honeycomb (bipartite): tr(H²) = 2 t² · n_NN_bonds = 2 t² · 3·Lx·Ly (Wallace 1947 / Castro Neto et al. 2009)",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end

@testset "Honeycomb TB — TightBindingChecksum vs real-space ED (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            lat = build_lattice(L2D_Honeycomb, Lx, Ly)
            H = build_tight_binding(lat, t)
            # A3 guard: if a future build_tight_binding change makes H non-Hermitian,
            # Symmetric(H) would silently use only the upper triangle and produce a
            # different effective matrix; fail loudly here instead.
            @assert norm(H - H') < 1e-12 "build_tight_binding returned non-Hermitian H at Lx=$(Lx), Ly=$(Ly), t=$(t)"
            ed_checksum = sum(abs2, eigvals(Symmetric(H)))
            verify(
                QAtlas.Honeycomb(),
                TightBindingChecksum(),
                Infinite();
                route=:ed_finite_size,
                independent=ed_checksum,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "ED black-box: Lattice2D.build_lattice(Honeycomb,Lx,Ly) + test/util/tight_binding.jl build_tight_binding → eigvals → Σ λ²",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end

@testset "Honeycomb TB — TightBindingMaxEnergy Γ-point literature pin (補完 after #449)" begin
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3), (3, 4))
        for t in (1.0, 1.5)
            verify(
                QAtlas.Honeycomb(),
                TightBindingMaxEnergy(),
                Infinite();
                route=:literature_value,
                independent=3.0 * t,
                agree_within=1e-10,
                at=["Lx=$(Lx)", "Ly=$(Ly)", "t=$(t)"],
                refs=[
                    "Honeycomb TB at Γ: |f(k=0)| = 3, so max(spectrum) = 3·|t| (independent of Lx, Ly for PBC; sensitive to chiral pair sign-flips that TightBindingChecksum cannot see)",
                ],
                fetch_kw=(; Lx=Lx, Ly=Ly, t=t),
            )
        end
    end
end
