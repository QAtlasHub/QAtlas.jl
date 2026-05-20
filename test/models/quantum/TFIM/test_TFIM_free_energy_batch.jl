# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/TFIM/test_TFIM_free_energy_batch.jl
#
# TFIM free-energy density f(β) at trivial limits.  Uses the existing
# thermal API with β = 1e6 to extract the T → 0 limit, where
# f → ε₀ (Pfeuty 1970 ground-state energy density).
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "TFIM — FreeEnergy trivial limits (#381 batch)" begin
    BETA = 1e6

    # /Infinite at three Pfeuty special points (matches PR #389 GSE values):
    #   h = 0 ⇒ f = -J ;  J = 0 ⇒ f = -h ;  h = J critical ⇒ f = -4J/π.
    for J in (0.5, 1.0, 2.0)
        verify(
            TFIM(; J=J, h=0.0),
            FreeEnergy(),
            Infinite();
            route=:second_closed_form,
            independent=-J,
            agree_within=1e-5,
            refs=["Pfeuty 1970: at h=0 dispersion flat Λ=2J ⇒ ε₀ = -J; β→∞ free energy → ε₀"],
            fetch_kw=(; beta=BETA),
        )
        verify(
            TFIM(; J=0.0, h=J),
            FreeEnergy(),
            Infinite();
            route=:second_closed_form,
            independent=-J,
            agree_within=1e-5,
            refs=["Pfeuty 1970: at J=0 dispersion flat Λ=2h ⇒ ε₀ = -h; β→∞ free energy → ε₀"],
            fetch_kw=(; beta=BETA),
        )
        verify(
            TFIM(; J=J, h=J),
            FreeEnergy(),
            Infinite();
            route=:second_closed_form,
            independent=-4 * J / π,
            agree_within=1e-5,
            refs=["Pfeuty 1970: at h=J critical point ε₀ = -4J/π; β→∞ free energy → ε₀"],
            fetch_kw=(; beta=BETA),
        )
    end

    # /OBC at J = 0: each spin contributes -h regardless of boundaries,
    # so f/N = -h exactly for any N.
    for h in (0.5, 1.0, 2.0)
        for N in (8, 12, 16)
            verify(
                TFIM(; J=0.0, h=h),
                FreeEnergy(),
                OBC(N);
                route=:second_closed_form,
                independent=-h,
                agree_within=1e-12,
                refs=["TFIM J=0 OBC: no Ising coupling ⇒ f/N = -h exactly for any N"],
                fetch_kw=(; beta=BETA),
            )
        end
    end

    # /OBC at h = 0: pure Ising chain with N-1 bonds on N sites ⇒
    # E₀/N = -(N-1)J/N exactly (boundary correction is the missing wrap-around bond).
    for J in (0.5, 1.0, 2.0)
        for N in (8, 12, 16)
            verify(
                TFIM(; J=J, h=0.0),
                FreeEnergy(),
                OBC(N);
                route=:second_closed_form,
                independent=-(N - 1) * J / N,
                agree_within=1e-5,
                refs=["TFIM h=0 OBC: N-1 bonds on N sites ⇒ ε₀ = -(N-1)J/N (boundary-correction closed form)"],
                fetch_kw=(; beta=BETA),
            )
        end
    end

    # /PBC at h = 0: pure Ising ring of N sites ⇒ ε₀ = -J exact (one bond per site).
    for J in (0.5, 1.0, 2.0)
        for N in (8, 12, 16)
            verify(
                TFIM(; J=J, h=0.0),
                FreeEnergy(),
                PBC(N);
                route=:second_closed_form,
                independent=-J,
                agree_within=1e-5,
                refs=["TFIM h=0 PBC: pure Ising ring, one bond per site ⇒ ε₀ = -J exactly"],
                fetch_kw=(; beta=BETA),
            )
        end
    end

    # /PBC at J = 0: same as OBC — f/N = -h exactly for any N.
    for h in (0.5, 1.0, 2.0)
        for N in (8, 12, 16)
            verify(
                TFIM(; J=0.0, h=h),
                FreeEnergy(),
                PBC(N);
                route=:second_closed_form,
                independent=-h,
                agree_within=1e-5,
                refs=["TFIM J=0 PBC: no Ising coupling ⇒ f/N = -h exactly"],
                fetch_kw=(; beta=BETA),
            )
        end
    end
end
