# =============================================================================
# TFIM Z-axis — verify() cards (WHY-correct plane)
#
# Split out of test/models/quantum/TFIM/test_TFIM_zaxis.jl (3.4 min on s05).
# Helpers _build_tfim_dense, _op_site, _SZ come from
# test/util/tfim_dense_ed.jl via runtests.jl ambient include.
# =============================================================================

@testset "TFIM z-axis — verification cards" begin
    # Pfeuty 1970 spontaneous magnetization: m_z = (1 - (h/J)²)^{1/8}
    # for h < J (ordered phase), and exactly 0 for h >= J.
    for (J, h) in ((1.0, 0.3), (1.0, 0.5), (1.0, 0.8))
        verify(
            TFIM(; J=J, h=h),
            MagnetizationZ(),
            Infinite();
            route=:second_closed_form,
            independent=(1 - (h / J)^2)^(1 / 8),
            agree_within=1e-9,
            refs=["Pfeuty 1970: m_z = (1 - (h/J)²)^{1/8} for h < J"],
        )
    end
    for (J, h) in ((1.0, 1.0), (1.0, 2.0))
        verify(
            TFIM(; J=J, h=h),
            MagnetizationZ(),
            Infinite();
            route=:limiting_case,
            independent=0.0,
            agree_within=1e-10,
            refs=["Pfeuty 1970: m_z = 0 for h >= J (disordered/critical)"],
        )
    end

    # ZZ static correlation at small N vs independent OBC dense ED
    let J = 1.0, h = 1.0, N = 6, i = 2, j = 5
        F = LinearAlgebra.eigen(_build_tfim_dense(N, J, h))
        ψ = F.vectors[:, 1]
        zz_ed = real(ψ' * (_op_site(_SZ, i, N) * (_op_site(_SZ, j, N) * ψ)))
        verify(
            TFIM(; J=J, h=h),
            ZZCorrelation(; mode=:static),
            OBC(N);
            route=:ed_finite_size,
            fetch_kw=(; i=i, j=j, beta=Inf),
            independent=zz_ed,
            agree_within=1e-8,
            refs=["Direct OBC dense-ED ⟨σz_i σz_j⟩ via _build_tfim_dense ground state"],
        )
    end
end
