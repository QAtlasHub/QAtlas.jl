# ─────────────────────────────────────────────────────────────────────────────
# test/models/classical/test_IsingSquare_PBC_thermo_ED_batch.jl
#
# ED-independent verify cards for IsingSquare/{Energy, SpecificHeat,
# ThermalEntropy}/PBC.  Direct brute-force enumeration over all 2^N spin
# configurations on the Lx × Ly torus (N=Lx·Ly), computing the canonical
# Z, E, F, S, C without touching the src transfer-matrix kernel.
#
# Capped at small (Lx, Ly) ∈ {(3,3)} for fast CI; L=2 dropped because the
# brute-force ED bond loop double-counts wrap-around bonds at L=2.
# Pure verify(); branches off main. Refs #381.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

function _brute_ising_thermo(Lx::Int, Ly::Int, J::Real, beta::Real)
    N = Lx * Ly
    Z = 0.0
    sumE = 0.0
    sumE2 = 0.0
    for s in 0:(2 ^ N - 1)
        spin = [((s >> k) & 1) == 1 ? 1.0 : -1.0 for k in 0:(N - 1)]
        E = 0.0
        for i in 1:Lx, j in 1:Ly
            idx = (i - 1) * Ly + j
            i2 = i % Lx + 1;
            idx2 = (i2 - 1) * Ly + j
            E += -J * spin[idx] * spin[idx2]
            j2 = j % Ly + 1;
            idx2 = (i - 1) * Ly + j2
            E += -J * spin[idx] * spin[idx2]
        end
        w = exp(-beta * E)
        Z += w
        sumE += E * w
        sumE2 += E * E * w
    end
    E_avg = sumE / Z
    F = -log(Z) / beta
    S = beta * (E_avg - F)
    C = beta^2 * (sumE2 / Z - E_avg^2)
    return E_avg / N, F / N, S / N, C / N
end

@testset "IsingSquare — {Energy, SpecificHeat, ThermalEntropy}/PBC vs brute-force ED (#381 batch)" begin
    # L=2 in any direction double-counts wrap-around bonds in the brute-force
    # ED loop (i % Lx + 1 collapses to the same neighbour), so we only sweep
    # sizes with L ≥ 3 in both directions.
    for (Lx, Ly) in ((3, 3),)
        for J in (0.5, 1.0)
            for beta in (0.3, 1.0)
                ed_E, ed_F, ed_S, ed_C = _brute_ising_thermo(Lx, Ly, J, beta)
                model = IsingSquare(; J=J, Lx=Lx, Ly=Ly)
                verify(
                    model,
                    Energy(:per_site),
                    PBC();
                    route=:ed_finite_size,
                    independent=ed_E,
                    at=["LxLy=$(Lx)x$(Ly)"],
                    agree_within=1e-7,
                    refs=[
                        "Brute-force ED: enumerate all 2^(Lx·Ly) configs, canonical ⟨H⟩_β/N on the torus",
                    ],
                    fetch_kw=(; beta=beta),
                )
                verify(
                    model,
                    SpecificHeat(),
                    PBC();
                    route=:ed_finite_size,
                    independent=ed_C,
                    at=["LxLy=$(Lx)x$(Ly)"],
                    agree_within=1e-7,
                    refs=["Brute-force ED: β²·(⟨E²⟩-⟨E⟩²)/N over all configs"],
                    fetch_kw=(; beta=beta),
                )
                verify(
                    model,
                    ThermalEntropy(),
                    PBC();
                    route=:ed_finite_size,
                    independent=ed_S,
                    at=["LxLy=$(Lx)x$(Ly)"],
                    agree_within=1e-7,
                    refs=["Brute-force ED: S/N = β·(E-F)/N over all configs"],
                    fetch_kw=(; beta=beta),
                )
            end
        end
    end
end
