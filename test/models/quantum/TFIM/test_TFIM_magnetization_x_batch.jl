# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_magnetization_x_batch.jl
#
# TFIM transverse magnetization <σ^x> at trivial limits:
#   * J = 0 (pure transverse field, ground state fully x-polarised) ⇒ m_x = sgn(h)
#   * h = 0 (pure Ising ferromagnet, ground states are z-eigenstates) ⇒ m_x = 0
# Verified across Infinite / OBC(N) / PBC(N) using the existing thermal API
# with beta → ∞ (here beta = 1e6). Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — MagnetizationX trivial limits (#381 batch)" begin
    BETA = 1e6  # β → ∞ ground-state extraction

    # J = 0, h > 0: pure transverse field, GS fully polarised ⇒ m_x = 1.
    for h in (0.5, 1.0, 2.0)
        # /Infinite
        verify(
            TFIM(; J=0.0, h=h),
            MagnetizationX(),
            Infinite();
            route=:second_closed_form,
            independent=1.0,
            agree_within=1e-12,
            refs=[
                "TFIM at J=0 is the pure transverse field; GS is the +x polarised product state ⇒ <σ^x> = 1",
            ],
            fetch_kw=(; beta=BETA),
        )
        # /OBC(N), /PBC(N) — same in any geometry at J=0 (no Ising coupling)
        for N in (8, 12, 16)
            verify(
                TFIM(; J=0.0, h=h),
                MagnetizationX(),
                OBC(N);
                route=:second_closed_form,
                independent=1.0,
                agree_within=1e-12,
                refs=[
                    "TFIM at J=0, OBC: BC-independent — GS still fully x-polarised ⇒ <σ^x> = 1",
                ],
                fetch_kw=(; beta=BETA),
            )
            verify(
                TFIM(; J=0.0, h=h),
                MagnetizationX(),
                PBC(N);
                route=:second_closed_form,
                independent=1.0,
                agree_within=1e-12,
                refs=[
                    "TFIM at J=0, PBC: BC-independent — GS still fully x-polarised ⇒ <σ^x> = 1",
                ],
                fetch_kw=(; beta=BETA),
            )
        end
    end

    # h = 0: pure Ising. With β finite the thermal trace is over Z₂-symmetric
    # states; the BdG-canonical-ensemble value of <σ^x> is identically zero
    # (no x-component in the Hamiltonian, no symmetry breaking from the
    # field channel).
    for J in (0.5, 1.0, 2.0)
        verify(
            TFIM(; J=J, h=0.0),
            MagnetizationX(),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=1e-12,
            refs=[
                "TFIM at h=0 (pure Ising): Hamiltonian has no x-component ⇒ <σ^x> = 0 (symmetry-unbroken thermal value)",
            ],
            fetch_kw=(; beta=BETA),
        )
        # /OBC(N), /PBC(N) at h = 0: finite-N TFIM GS preserves Z₂ exactly,
        # so <σ^x> = 0 by parity for every N (no spontaneous breaking at
        # finite N). Sweep N consistent with the J=0 branch above.
        for N in (8, 12, 16)
            verify(
                TFIM(; J=J, h=0.0),
                MagnetizationX(),
                OBC(N);
                route=:second_closed_form,
                independent=0.0,
                agree_within=1e-12,
                refs=["TFIM at h=0, OBC: exact Z₂ symmetry at finite N ⇒ <σ^x> = 0"],
                fetch_kw=(; beta=BETA),
            )
            verify(
                TFIM(; J=J, h=0.0),
                MagnetizationX(),
                PBC(N);
                route=:second_closed_form,
                independent=0.0,
                agree_within=1e-12,
                refs=["TFIM at h=0, PBC: exact Z₂ symmetry at finite N ⇒ <σ^x> = 0"],
                fetch_kw=(; beta=BETA),
            )
        end
    end
end
