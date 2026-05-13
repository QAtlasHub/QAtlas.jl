# ─────────────────────────────────────────────────────────────────────────────
# BCFT — Boundary Conformal Field Theory with Cardy boundary states.
#
# For a 2-D CFT defined on a strip (or upper half plane) with conformal
# boundary conditions, Cardy (1989) classified the consistent boundary
# states `|a⟩` in one-to-one correspondence with the chiral primaries
# of the bulk CFT.  Each Cardy state has an associated Affleck-Ludwig
# `g`-function (Affleck-Ludwig 1991):
#
#     g_a = ⟨0|a⟩ = S_{a,0} / √(S_{0,0}),
#
# where `S_{a,b}` is the modular S-matrix of the bulk CFT.  The
# **boundary entropy** is then
#
#     s_bdy = log g_a,
#
# which is non-decreasing along boundary RG flows — the g-theorem,
# proven by Friedan-Konechny (2004).
#
# Phase 1 ships only the **critical Ising CFT** (M(4,3), c = 1/2) Cardy
# boundary states, labelled by the three Ising primaries:
#
# | primary  | h    | Cardy state            | g          |
# | -------- | ---- | ---------------------- | ---------- |
# | 1        | 0    | "fixed +" / "fixed −"  | 1/√2       |
# | ε        | 1/2  | "free"                 | 1          |
# | σ        | 1/16 | "Cardy superposition"  | 1/√2       |
#
# Other CFTs (minimal models, free boson, WZW, …) and richer boundary
# observables (boundary entropy along RG flows, defect g-functions)
# are tracked as **Phase 2**.
#
# References:
#   - J. L. Cardy, Nucl. Phys. B 324, 581 (1989).
#   - I. Affleck, A. W. W. Ludwig, Phys. Rev. Lett. 67, 161 (1991).
#   - D. Friedan, A. Konechny, Phys. Rev. Lett. 93, 030402 (2004).
# ─────────────────────────────────────────────────────────────────────────────

"""
    BCFT() <: AbstractQAtlasModel

Boundary Conformal Field Theory with Cardy (1989) boundary states and
Affleck-Ludwig (1991) `g`-functions.

**Phase 1 ships only Ising-CFT Cardy boundary states**; other CFTs
(minimal models, free boson, WZW, …) are tracked as Phase 2.

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                                       |
| ------------------------------ | ---------- | -------------------------------------------- |
| [`ResidualEntropy`](@ref)      | `Infinite` | analytic `log g` for Ising Cardy states      |

The Cardy state is selected via the `state` keyword to `fetch`:

- `state=:fixed`, `:fixed_plus`, `:fixed_minus` ⟹ `g = 1/√2`, `log g = -log(2)/2`
- `state=:free`                                  ⟹ `g = 1`,    `log g = 0`
- `state=:sigma`                                 ⟹ `g = 1/√2`, `log g = -log(2)/2`

# References

- J. L. Cardy, *Nucl. Phys. B* **324**, 581 (1989).
- I. Affleck, A. W. W. Ludwig, *Phys. Rev. Lett.* **67**, 161 (1991).
- D. Friedan, A. Konechny, *Phys. Rev. Lett.* **93**, 030402 (2004).
"""
struct BCFT <: AbstractQAtlasModel end

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 1: Ising Cardy boundary entropy log g via ResidualEntropy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::BCFT, ::ResidualEntropy, ::Infinite; state::Symbol=:fixed) -> Float64

Return the Affleck-Ludwig boundary entropy `log g_a` for an Ising-CFT
Cardy boundary state `a`.

| `state`                                  | g       | `log g`        |
| ---------------------------------------- | ------- | -------------- |
| `:fixed`, `:fixed_plus`, `:fixed_minus`  | `1/√2`  | `-log(2)/2`    |
| `:free`                                  | `1`     | `0`            |
| `:sigma`                                 | `1/√2`  | `-log(2)/2`    |

Throws a `DomainError` for any other `state` symbol (Phase 1 only
supports Ising Cardy boundary states; other CFTs are Phase 2).
"""
function fetch(::BCFT, ::ResidualEntropy, ::Infinite; state::Symbol=:fixed, kwargs...)
    if state === :fixed || state === :fixed_plus || state === :fixed_minus
        return -log(2) / 2   # log(1/√2) = -(1/2) log 2
    elseif state === :free
        return 0.0            # log 1
    elseif state === :sigma
        return -log(2) / 2   # log(1/√2)
    else
        throw(DomainError(state,
            "BCFT ResidualEntropy: Phase 1 supports Ising Cardy states only " *
            "(:fixed_plus, :fixed_minus, :free, :sigma). Got state=:$state. " *
            "Other CFTs and boundary states tracked as Phase 2."))
    end
end
