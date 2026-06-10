# Verification — putting the physics in correctly

A value in `src/` is rigorous only once an **independent** route reproduces it.
This is the heart of QAtlas: the stored number and its check must not share a
derivation.

## The verify card

Every new `(model, quantity, bc)` triple needs at least one `verify(...)` card
in `test/models/<class>/<Model>/...`. The signature is frozen in
`test/util/verify.jl`:

```julia
verify(model, quantity, bc;
       route::Symbol,            # one of the allowed routes
       independent,             # scalar OR convergence vector (the reference value)
       agree_within::Real,      # ABSOLUTE tolerance
       refs::Vector{<:AbstractString},
       fetch_kw::NamedTuple = (;),   # parameters passed into the black-box fetch
       at = nothing,            # x-axis for convergence vectors
       subject_extract = nothing)
```

Two properties make the card non-tautological:

- The subject is **fetched inside** `verify` — you cannot precompute it and pass
  it in.
- The `independent` value must be derivable **without** running `fetch`'s code
  path. Re-running the same `quadgk` integral is **not** independent.

## Choosing an independent route

| Route | Independent because… |
|-------|----------------------|
| `:ed_finite_size` | exact diagonalisation at small `N` (optionally `N→∞` extrapolated) — a completely different computation from a closed form |
| `:second_closed_form` | a different analytic derivation of the same quantity |
| `:limiting_case` | a known value at a special point (`T=0`, `Δ=0`, `β→0/∞`) |
| `:delegation_invariant` | model X at parameter p *is exactly* model Y (e.g. `XXZ1D(Δ=0)` ≡ free fermion) |
| `:literature_value` | a published DMRG/MC number (cite the DOI + table) |
| `:sum_rule` | an independent analytic identity |

Pick the strongest route available. Limiting cases alone are weak; prefer ED or
a second closed form when they exist.

## The strongest checks

- **Independent Boltzmann-ED.** For a thermal quantity computed from cumulants /
  a closed form, diagonalise a finite ring and sum `Z = Σ e^{-βE}` directly.
  This shares the Hamiltonian but nothing else, so it catches every bug in the
  closed-form *assembly* (signs, factors, the cumulant→thermo algebra).

- **Cross-model exact match.** When two models coincide at a dual/limit point,
  check one model's new quantity against the *other model's independent
  implementation*. Example: at the Jordan–Wigner-dual point
  `Kitaev1D(μ=−2h, t=J, Δ=J)` has `E(k) ≡ Λ_k` of `TFIM(J,h)`, so the Kitaev
  free energy must equal the TFIM free energy to ~1e-6 — this even pins
  zero-point / ½-doubling constants that limits and autodiff cannot.

- **Literature cross-check via `doiget`.** For literature-grounded coefficients,
  the reference value must come from the actual paper, not your derivation. See
  [citations.md](citations.md).

## What internal consistency is *not*

Autodiff identities (`c_v = β²∂²lnZ`, `s = β(ε−f)`), scaling relations
(Rushbrooke, Widom), and "two of my own routes agree" are **necessary but not
sufficient**. They confirm the algebra is self-consistent; they cannot detect a
convention error baked into both sides. Always pair them with at least one
externally-anchored route (ED, cross-model, or literature).

## Tests must hit the real code path

Do not ship a primitive validated only by a toy 5×5 unit test. Exercise the
real caller with a dense / independent reference (an integration test), not just
the helper in isolation. A green single-file run can also hide collisions — run
the multi-file suite before pushing (a util method that is *more specific* than
an existing one silently hijacks dispatch).

## Plain tests still have a place

Exception-shape checks (`@test_throws DomainError` on `β ≤ 0`) and
identity-between-fetched-values checks (Gibbs `s = β(u−f)`) don't fit the
per-quantity card schema — keep them as plain `@test` / `@test_throws`, as a
complement to (not a replacement for) verify cards.
