# ─────────────────────────────────────────────────────────────────────────────
# test/universalities/test_universality_central_charge_lit.jl
#
# Literature-value pins for Universality{X}/CentralCharge/Infinite.
#
# This file restores a piece of the WHY-plane coverage removed in
# PR #449 (zero-legacy phase 1): the deleted
# test/verification/universality/test_universality_literature_values.jl
# and test/verification/universality/test_entanglement_central_charge.jl
# both pinned the central charges of the five Universality{C} classes
# against the literature decimals cited in src/universalities/. The
# new file does the same job through the verify() framework — each
# card is a structural literature_value cross-check that catches any
# accidental drift of a stored c constant.
#
# Hubs registered (all new — INVENTORY had zero Universality/CentralCharge
# entries before this PR):
#   Universality(:Ising)/CentralCharge/Infinite       (BPZ 1984)
#   Universality(:Potts3)/CentralCharge/Infinite      (Dotsenko 1984)
#   Universality(:Potts4)/CentralCharge/Infinite      (di Francesco-Mathieu-Sénéchal §12.3)
#   Universality(:XY)/CentralCharge/Infinite          (Kosterlitz 1974)
#   Universality(:Heisenberg)/CentralCharge/Infinite  (Affleck-Haldane 1987)
#
# Pure verify(); branches off main. Refs #381 / restores coverage lost in #449.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality CentralCharge — literature pins (補完 after #449)" begin
    verify(
        Universality(:Ising),
        CentralCharge(),
        Infinite();
        route=:literature_value,
        independent=1//2,
        agree_within=0,
        refs=[
            "Belavin-Polyakov-Zamolodchikov 1984 (Nucl. Phys. B 241, 333): 2D Ising CFT = M(4,3), c = 1/2",
        ],
        fetch_kw=(; d=2),
    )

    verify(
        Universality(:Potts3),
        CentralCharge(),
        Infinite();
        route=:literature_value,
        independent=4//5,
        agree_within=0,
        refs=[
            "Dotsenko 1984 (Nucl. Phys. B 235, 54); di Francesco-Mathieu-Sénéchal §7.4: 2D 3-state Potts = M(6,5), c = 4/5",
        ],
        fetch_kw=(; d=2),
    )

    verify(
        Universality(:Potts4),
        CentralCharge(),
        Infinite();
        route=:literature_value,
        independent=1//1,
        agree_within=0,
        refs=[
            "di Francesco-Mathieu-Sénéchal 1997 (Conformal Field Theory) §12.3: 2D 4-state Potts at marginal compact-boson self-dual radius, c = 1",
        ],
        fetch_kw=(; d=2),
    )

    verify(
        Universality(:XY),
        CentralCharge(),
        Infinite();
        route=:literature_value,
        independent=1//1,
        agree_within=0,
        refs=[
            "Kosterlitz 1974 (J. Phys. C 7, 1046); di Francesco-Mathieu-Sénéchal §6: 2D XY (BKT critical line) = free compact boson, c = 1",
        ],
        fetch_kw=(; d=2),
    )

    verify(
        Universality(:Heisenberg),
        CentralCharge(),
        Infinite();
        route=:literature_value,
        independent=1//1,
        agree_within=0,
        refs=[
            "Affleck-Haldane 1987 (Phys. Rev. B 36, 5291): spin-1/2 Heisenberg chain in SU(2)_1 WZW universality, c = 1",
        ],
        fetch_kw=(; d=1),
    )
end
