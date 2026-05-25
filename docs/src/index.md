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

All reports go to [sotashimozono/QAtlas.jl Issues](https://github.com/sotashimozono/QAtlas.jl/issues).
Corrections and pull requests are equally welcome.
