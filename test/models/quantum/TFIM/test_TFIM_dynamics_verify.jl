# =============================================================================
# TFIM dynamics — verify() cards (WHY-correct plane)
#
# Split out of the legacy test_TFIM_dynamics.jl (PR vs next: refactor for
# CI shard balance — original file was ~17 min on s14 because Layers 2b+2c
# ran heavy BdG Pfaffian sweeps at N=200-240 inside the same file).
#
# Helpers _build_tfim_dense, _op_site, _SZ, _SX come from test/util/
# tfim_dense_ed.jl via runtests.jl ambient include.
# =============================================================================

@testset "TFIM dynamics — verification cards" begin
    # t = 0: the dynamic ZZ correlator reduces to the static ground-state
    # correlator (independent route: direct dense ED of the GS).
    let J = 1.0, h = 1.0, N = 6, i = 2, j = 4
        F = eigen(_build_tfim_dense(N, J, h))
        ψ = F.vectors[:, 1]
        zz0 = real(ψ' * (_op_site(_SZ, i, N) * (_op_site(_SZ, j, N) * ψ)))
        verify(
            TFIM(; J=J, h=h),
            ZZCorrelation(; mode=:dynamic),
            OBC(N);
            route=:limiting_case,
            fetch_kw=(; i=i, j=j, t=0.0),
            independent=zz0,
            agree_within=1e-8,
            refs=["t=0 limit: ⟨σz_i(0) σz_j⟩ = static GS correlator (dense ED)"],
        )
    end

    # GS energy is time-independent (conserved): equals min eigenvalue
    let J = 1.0, h = 1.3, N = 6
        verify(
            TFIM(; J=J, h=h),
            Energy(),
            OBC(N);
            route=:ed_finite_size,
            fetch_kw=(; beta=Inf),
            independent=dense_spectrum(_build_tfim_dense(N, J, h))[1],
            agree_within=1e-9,
            refs=["GS energy = min eigenvalue of _build_tfim_dense"],
        )
    end
end
