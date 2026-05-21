# =============================================================================
# KitaevHoneycomb - verify() cards + MassGap closed form (testsets 11-12)
#
# Split out of test/models/quantum/KitaevHoneycomb/test_KitaevHoneycomb.jl
# (2.5 min on s06) so each top-level group runs on its own shard.
# =============================================================================

using QAtlas: KitaevHoneycomb, Energy, MassGap, Infinite, PBC, OBC
using Test
using LinearAlgebra: eigvals, Hermitian, I as I_mat
using SparseArrays: sparse, spzeros
using KrylovKit: eigsolve
using Lattice2D: build_lattice, Honeycomb, OpenAxis, PeriodicAxis
using LatticeCore: num_sites, bonds

# Reference value taken from Baskaran, Mandal & Shankar (PRL 98, 247201)
# and reproduced in Feng et al. DMRG benchmarks of the Kitaev model:
# at the isotropic point (Kx = Ky = Kz = 1) the ground state energy
# per site for H = -Σ Kγ σγ σγ is
#
#     ε₀ ≈ -0.7872986...
#
# i.e. `-⟨|f(k)|⟩ / 2` in the conventions of the module.

const ε_isotropic_TL = -0.7872986216706852

@testset "KitaevHoneycomb — verification cards" begin
    # Isotropic Kitaev honeycomb ground-state energy density: the exact
    # Kitaev-2006 Brillouin-zone integral value (literature constant).
    verify(
        KitaevHoneycomb(; Kx=1.0, Ky=1.0, Kz=1.0),
        Energy(),
        Infinite();
        route=:literature_value,
        independent=-0.7872986216706852,
        agree_within=1e-6,
        refs=["Kitaev 2006: isotropic honeycomb e0 ≈ -0.78729862 |K|"],
    )
end

# ── additional verification card (#381 batch) ─────────────────────────────
@testset "KitaevHoneycomb — MassGap/Infinite closed form (#381 batch)" begin
    # Kitaev 2006: Δ = 2 · max(|K_max| − |K_others_sum|, 0).
    # In gapless A/B/C phase (|K_γ| ≤ sum of others for all γ) ⇒ Δ = 0.
    # In gapped A_γ phase (one |K_γ| exceeds the sum of the other two) the
    # excess sets the single-Majorana gap.
    for (Kx, Ky, Kz, Δ_expected) in (
        # Gapless isotropic / B-phase points: Δ = 0
        (1.0, 1.0, 1.0, 0.0),
        (1.0, 1.0, 1.5, 0.0),     # gapless interior: |Kz|=1.5 < 1+1 = 2
        (0.5, 0.5, 0.9, 0.0),     # gapless interior: |Kz|=0.9 < 0.5+0.5 = 1.0
        # True B-phase boundary: |K_max| = sum_of_others (Δ = 0 by the max(·,0) clamp).
        (1.0, 1.0, 2.0, 0.0),     # boundary: |Kz|=2.0 = 1+1
        (0.5, 0.5, 1.0, 0.0),     # boundary: |Kz|=1.0 = 0.5+0.5
        # Gapped A_z phase: Kz > Kx + Ky
        (1.0, 1.0, 3.0, 2 * (3.0 - 2.0)),  # = 2.0
        (0.5, 0.5, 2.0, 2 * (2.0 - 1.0)),  # = 2.0
        # Gapped A_x phase: Kx > Ky + Kz
        (3.0, 1.0, 1.0, 2.0),
        # Gapped A_y phase: Ky > Kx + Kz
        (0.5, 4.0, 0.5, 2 * (4.0 - 1.0)),  # = 6.0
    )
        verify(
            KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz),
            MassGap(),
            Infinite();
            route=:second_closed_form,
            independent=Δ_expected,
            agree_within=1e-12,
            refs=[
                "Kitaev 2006 Annals 321: Δ = 2·max(|K_max| − |K_others_sum|, 0); gapped iff one |K_γ| > sum of others",
            ],
        )
    end
end
