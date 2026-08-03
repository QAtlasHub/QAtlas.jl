# QAtlas.jl

**QAtlas** (QUAntum Reference Table for Exact Tests) is a curated
dictionary of rigorous results in quantum and statistical physics.
Every stored value is traced to a specific publication and
cross-validated against independent calculations.

## Key Features

Unlike typical numerical libraries, QAtlas focuses on **authoritative
reference values** — exact analytical results, high-precision
conformal bootstrap bounds, and Bethe ansatz solutions. Each value
is accompanied by:

1. **Precise citation**: author, year, journal, equation number
2. **Derivation sketch**: enough to independently verify
3. **Cross-validation**: tested against independent computation
4. **Connections**: linked to universality classes and other models

### Two kinds of stored value, and the registry says which

The paragraph above describes most of the atlas but **not all of it**, and the
difference is worth stating rather than leaving a reader to discover.

A row registered `cost = :exponential` is many-body **dense exact
diagonalisation** on a `d^N` Hilbert space, with a declared `max_size` (12 for
spin-½).  Measured on the loaded registry: **69 of 518 rows (13.3%)**,
concentrated in three models — `S1Heisenberg1D` (22 of 24), `Heisenberg1D` (23 of
35), `XXZ1D` (23 of 37), where they are the whole `OBC` surface.

Those rows are **not** "authoritative reference values" in the sense above.
There is no publication for the local magnetisation of the eight-site open
Heisenberg chain at `β = 1`, and any consumer with an ED routine reproduces the
number in minutes.  What they are is a **cross-validation instrument**:

- the small-`N` ground truth an analytic closed form is *checked against* (the
  TFIM entanglement entropy documents exactly this — "matches the full-ED
  reference at every small `N`, verified to 1e-10");
- the target a new DMRG / MPS / QMC implementation validates itself on before
  it is trusted at sizes where nothing exact exists.

Two things that are easy to assume and are false, both measured:

- **They do not duplicate the closed forms.** Every `(model, quantity)` pair
  carrying both has the ED row at a finite `OBC` chain and the analytic row at
  `Infinite` — the thermodynamic limit is a *different question*, not a second
  route to the same number.  No `(model, quantity, bc)` triple carries both, and
  a test asserts it stays that way.  For **56** `(model, quantity)` pairs, ED is
  the only route there is.
- **The cost is not a detail.** `cost = :exponential` requires `max_size` at
  registration precisely so a route cannot hide where it stops, and the limit
  bites in practice: one region-entropy bag costs 0.0 s on `TFIM/OBC` and
  **338.7 s** on `S1Heisenberg1D/OBC`, which is why the latter carries an
  exclusion in `region_registry.jl`.

So: read `cost` alongside `status`.  `status = :exact` says the value is right;
`cost = :exponential` says it is an instrument you could have built yourself,
and tells you the size beyond which it does not exist.

## Quick Start

```julia
using QAtlas

# Onsager critical temperature
Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature())

# TFIM ground-state energy
E₀ = QAtlas.fetch(:TFIM, :energy, OBC(); N=16, J=1.0, h=0.5)

# 2D Ising universality: exact exponents (Rational)
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
# (β = 1//8, ν = 1//1, γ = 7//4, η = 1//4, ...)
```

## Contents

- **[Models](models/index.md)** — exact solutions for classical and quantum models
- **[Universality Classes](universalities/index.md)** — critical exponents and scaling relations across dimensions
- **[Verification](verification/index.md)** — five-layer testing strategy ensuring physical correctness
- **[Methods](methods/index.md)** — computational techniques with physical justification
- **[Derivation Notes](calc/jw-tfim-bdg.md)** — step-by-step calculations
- **[API Reference](api.md)** — full Julia docstring index

## Reporting Errors

Every page has a **Report an issue** button fixed at the top-right of the
screen. Clicking it opens a pre-filled GitHub issue with the current page
URL — no copy-paste needed.

Individual sections also show a small **report** link when you hover over
an H2 or H3 heading. Use it to flag a specific derivation or formula that
looks wrong.

All reports go to [QAtlasHub/QAtlas.jl Issues](https://github.com/QAtlasHub/QAtlas.jl/issues).
Corrections and pull requests are equally welcome.
