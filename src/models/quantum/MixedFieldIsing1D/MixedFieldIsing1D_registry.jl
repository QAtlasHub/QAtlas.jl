# models/quantum/MixedFieldIsing1D/MixedFieldIsing1D_registry.jl
#
# One `@register` line per natively-implemented (model, quantity, bc)
# triple for the MixedFieldIsing1D model.  See `src/core/registry.jl`
# for the metadata schema and `TFIM_registry.jl` for the canonical
# pattern.

@register(
    MixedFieldIsing1D,
    MassGap,
    Infinite,
    method=:delegation,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_mixed_field_ising1d.jl",
    references=["Pfeuty 1970", "McCoy-Wu 1978"],
    notes="h_z = 0 delegated to TFIM (closed-form Δ = 2|h_x − J|); h_z ≠ 0 non-integrable, deferred to Phase 2 (ED/DMRG).",
)
