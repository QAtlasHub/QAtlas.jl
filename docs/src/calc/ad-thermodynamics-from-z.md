# Thermodynamic Quantities from ``\partial \ln Z`` via Automatic Differentiation

## Main result

For a classical or quantum statistical-mechanical system with
partition function ``Z(\beta, \{g_{i}\})`` — a smooth function of
the inverse temperature ``\beta = 1/T`` and couplings
``\{g_{i}\}`` (e.g. ``J``, ``h``) — all equilibrium thermodynamic
observables are derivatives of ``\ln Z``:

```math
\boxed{\;
\begin{aligned}
F(\beta) \;&=\; -\,\tfrac{1}{\beta}\ln Z, & &\text{(free energy)}\\[2pt]
\langle E\rangle \;&=\; -\,\frac{\partial\ln Z}{\partial\beta}
                   \;=\; \frac{\partial(\beta F)}{\partial\beta}, & &\text{(internal energy)}\\[2pt]
\langle E^{2}\rangle - \langle E\rangle^{2}
 \;&=\; \frac{\partial^{2}\ln Z}{\partial\beta^{2}}
 \;=\; -\,\frac{\partial\langle E\rangle}{\partial\beta}, & &\text{(energy variance)}\\[2pt]
C_{v}(\beta) \;&=\; k_{B}\beta^{2}\bigl(\langle E^{2}\rangle
                                       - \langle E\rangle^{2}\bigr)
                    \;=\; k_{B}\beta^{2}\frac{\partial^{2}\ln Z}{\partial\beta^{2}}, & &\text{(specific heat)}\\[2pt]
S(\beta) \;&=\; -\,\frac{\partial F}{\partial T}
               \;=\; k_{B}\bigl(\ln Z + \beta\langle E\rangle\bigr), & &\text{(entropy)}\\[2pt]
\langle g_{i}\text{-conjugate}\rangle
 \;&=\; -\,\tfrac{1}{\beta}\frac{\partial\ln Z}{\partial g_{i}}. & &\text{(coupling derivative)}
\end{aligned}
\;}
```

(Setting ``k_{B} = 1`` throughout for notational brevity.)

All seven formulas are elementary consequences of the
Boltzmann-weight identity $\partial_{\beta} e^{-\beta E_{s}} =
-E_{s}\,e^{-\beta E_{s}}`` plus the chain rule. When ``Z$ is
implemented numerically — e.g. as ``Z = \mathrm{Tr}(T^{L_{x}})`` for
a transfer matrix or as ``Z = \sum_{\{s\}}e^{-\beta H(\{s\})}`` for a
brute-force enumeration — all of these observables become
**automatically differentiable** via dual-number arithmetic
(ForwardDiff.jl), provided the numerical pipeline itself composes
only algebraic and analytic operations (matrix multiplication,
`log`, `tr`, `sum`, etc.).

In particular, QAtlas's `IsingSquare` transfer-matrix path uses
``Z = \mathrm{Tr}(T^{L_{x}})`` via repeated squaring — **not** an
eigenvalue decomposition — because LAPACK does not accept
dual-number inputs. All IsingSquare thermodynamic APIs
(`FreeEnergy`, `SpecificHeat`, `ThermalEntropy`, …) dispatch
through `ForwardDiff.derivative` on ``\ln Z(\beta)`` and are verified
against brute-force ensemble averages to ``\sim 10^{-14}`` relative
error (see `test/verification/test_ising_ad_thermodynamics.jl`).

---

## Setup

### Statistical-mechanical background

Given a Hamiltonian ``H(\{s\})`` over a state space ``\{s\}``, the
canonical partition function at inverse temperature ``\beta`` is

```math
Z(\beta) \;=\; \sum_{\{s\}} e^{-\beta H(\{s\})}.
```

The probability of a state is ``p_{s}(\beta) = e^{-\beta H_{s}}/Z``.
Thermal averages are ``\langle A\rangle = \sum_{s} A(s)\,p_{s}``. The
free energy is ``F = -T\ln Z = -\tfrac{1}{\beta}\ln Z``.

### Implementation target

For the 2D classical Ising model on an ``L_{x}\times L_{y}`` square
lattice with PBC and coupling ``J``, QAtlas computes

```math
Z_{\rm Ising}(\beta, J; L_{x}, L_{y})
 \;=\; \mathrm{Tr}\bigl(T^{L_{x}}\bigr),
```

via the symmetric transfer matrix ``T`` derived in
[`transfer-matrix-symmetric-split`](transfer-matrix-symmetric-split.md).
The goal is to obtain $F, \langle E\rangle, C_{v}, S,
\langle\sigma\sigma\rangle$ as (higher-order) derivatives of
``\ln Z`` using ForwardDiff.

### Goal

(i) Derive every Main-result identity from the definition of ``Z``
and elementary calculus.
(ii) Explain why forward-mode automatic differentiation via dual
numbers carries these identities through the numerical
implementation.
(iii) Specify the implementation constraints
(``Z = \mathrm{Tr}(T^{L_{x}})`` via matmul + `tr`, no LAPACK
eigensolvers).

---

## Calculation

### Step 1 — ``\langle E\rangle = -\partial\ln Z/\partial\beta``

Differentiate ``Z = \sum_{s} e^{-\beta H_{s}}`` with respect to
``\beta``:

```math
\frac{\partial Z}{\partial\beta}
 \;=\; -\sum_{s} H_{s}\,e^{-\beta H_{s}}
 \;=\; -Z\,\langle H\rangle,
```

using ``\langle H\rangle = Z^{-1}\sum_{s} H_{s}\,e^{-\beta H_{s}}``.
Divide by ``Z``:

```math
\frac{1}{Z}\frac{\partial Z}{\partial\beta}
 \;=\; \frac{\partial\ln Z}{\partial\beta}
 \;=\; -\,\langle E\rangle.
\tag{1}
```

Equivalently ``\langle E\rangle = -\partial\ln Z/\partial\beta``
(identifying ``E = H`` at thermal equilibrium).

### Step 2 — Variance identity $\partial^{2}\ln Z/\partial\beta^{2}
             = \langle E^{2}\rangle - \langle E\rangle^{2}$

Differentiate (1) once more:

```math
\frac{\partial^{2}\ln Z}{\partial\beta^{2}}
 \;=\; -\,\frac{\partial\langle E\rangle}{\partial\beta}
 \;=\; -\,\frac{\partial}{\partial\beta}
    \!\left[\frac{\sum_{s} H_{s}\,e^{-\beta H_{s}}}{Z}\right].
```

Apply the quotient rule:

```math
\frac{\partial}{\partial\beta}
  \!\left[\frac{\sum_{s} H_{s}\,e^{-\beta H_{s}}}{Z}\right]
 \;=\; \frac{1}{Z}\,\frac{\partial}{\partial\beta}
       \!\sum_{s} H_{s}\,e^{-\beta H_{s}}
    \;-\; \frac{1}{Z^{2}}\,
        \bigl(\sum_{s} H_{s}\,e^{-\beta H_{s}}\bigr)
        \,\frac{\partial Z}{\partial\beta}.
```

The numerator of the first term is $\sum_{s}(-H_{s}^{2})\,
e^{-\beta H_{s}} = -Z\,\langle H^{2}\rangle$. The second term is
$-\langle E\rangle \cdot \langle E\rangle\cdot(-Z)/Z = \langle
E\rangle^{2}`` (using ``\partial Z/\partial\beta = -Z\langle E\rangle$).

Putting it together,

```math
\frac{\partial\langle E\rangle}{\partial\beta}
 \;=\; -\langle E^{2}\rangle + \langle E\rangle^{2},
```

and therefore

```math
\boxed{\;
\frac{\partial^{2}\ln Z}{\partial\beta^{2}}
 \;=\; \langle E^{2}\rangle - \langle E\rangle^{2}
 \;=\; \mathrm{Var}(E) \;\ge\; 0.
\;}
\tag{2}
```

This is the **fluctuation-dissipation identity** in the
thermodynamic form: the variance of ``E`` equals the second
logarithmic derivative of ``Z`` with respect to ``\beta``. Positivity
``\mathrm{Var}(E) \ge 0`` corresponds to the convexity of ``-\ln Z``
in ``\beta`` — a thermodynamic stability requirement.

### Step 3 — Specific heat ``C_{v} = \beta^{2}\mathrm{Var}(E)``

The specific heat at constant volume is defined by

```math
C_{v} \;\equiv\; \frac{\partial\langle E\rangle}{\partial T}
 \;=\; \frac{d\beta}{dT}\,\frac{\partial\langle E\rangle}{\partial\beta}
 \;=\; -\frac{1}{T^{2}}\,\frac{\partial\langle E\rangle}{\partial\beta}
 \;=\; -\,k_{B}\beta^{2}\,\frac{\partial\langle E\rangle}{\partial\beta},
```

using ``\beta = 1/(k_{B} T)``. Substituting (2):

```math
\boxed{\;
C_{v} \;=\; k_{B}\beta^{2}\,\mathrm{Var}(E)
        \;=\; k_{B}\beta^{2}\,\frac{\partial^{2}\ln Z}{\partial\beta^{2}}.
\;}
\tag{3}
```

### Step 4 — Free energy and entropy

The Helmholtz free energy is ``F = -\beta^{-1}\ln Z``, so

```math
\frac{\partial(\beta F)}{\partial\beta}
 \;=\; -\,\frac{\partial\ln Z}{\partial\beta}
 \;=\; \langle E\rangle,
```

consistent with (1). The entropy is ``S = -\partial F/\partial T``.
Using ``F = -T\ln Z`` and the chain rule $\partial/\partial T =
-\beta^{2}\,\partial/\partial\beta\,\cdot(1/k_{B})`` (with ``k_{B} =
1$ henceforth),

```math
S \;=\; -\frac{\partial F}{\partial T}
 \;=\; \frac{\partial(T\ln Z)}{\partial T}
 \;=\; \ln Z + T\,\frac{\partial\ln Z}{\partial T}
 \;=\; \ln Z - T\beta^{2}\frac{\partial\ln Z}{\partial\beta}
 \;=\; \ln Z + \beta\langle E\rangle.
```

(Using ``T\beta^{2} = \beta`` and (1) at the last step.) This is the
standard Gibbs relation ``F = \langle E\rangle - T S``,
rearranged.

### Step 5 — Coupling-derivative identities

For any coupling ``g_{i}`` entering ``H`` through
``H(\{s\}) = H_{0}(\{s\}) + g_{i}\,\mathcal{O}_{i}(\{s\})``,

```math
\frac{\partial Z}{\partial g_{i}}
 \;=\; -\beta\sum_{s}\mathcal{O}_{i}(s)\,e^{-\beta H_{s}}
 \;=\; -\beta\,Z\,\langle\mathcal{O}_{i}\rangle,
```

so

```math
\boxed{\;\langle\mathcal{O}_{i}\rangle
 \;=\; -\,\tfrac{1}{\beta}\,\frac{\partial\ln Z}{\partial g_{i}}.\;}
\tag{4}
```

In the Ising case, setting ``g = J`` with $\mathcal{O}_{J} = -\sum
s_{u}s_{v}$ gives
``\langle\sum s_{u}s_{v}\rangle = \beta^{-1}\partial_{J}\ln Z``, i.e.
the nearest-neighbour correlator $\langle s_{u}s_{v}\rangle =
\beta^{-1}\partial_{J}\ln Z / N_{\rm bonds}``. Setting ``g = h$ with
``\mathcal{O}_{h} = -\sum s_{i}`` gives the magnetisation
``\langle M\rangle = \beta^{-1}\partial_{h}\ln Z``.

### Step 6 — Automatic differentiation via dual numbers

Forward-mode AD implements ``\partial f/\partial x`` by replacing
real inputs with **dual numbers** ``x = a + b\,\varepsilon`` where
``\varepsilon^{2} = 0`` (a formal infinitesimal). Any elementary
operation then propagates the dual part linearly:

```math
f(a + b\varepsilon) \;=\; f(a) + b\,f'(a)\,\varepsilon + O(\varepsilon^{2})
 \;=\; f(a) + b\,f'(a)\,\varepsilon,
```

using ``\varepsilon^{2} = 0``. Starting with $\beta = \beta_{0} +
1\cdot\varepsilon`` and computing ``\ln Z(\beta)$ via the full
numerical pipeline, the dual part of the result equals
``\partial\ln Z/\partial\beta|_{\beta_{0}}`` exactly — no
finite-difference approximation, no round-off in the derivative.

For **second** derivatives, iterate: use a dual-of-dual (or
hyperdual) type. ForwardDiff.jl implements this via

```julia
ForwardDiff.derivative(b -> ForwardDiff.derivative(b2 -> log(Z(b2)), b), β)
```

which returns ``\partial^{2}\ln Z/\partial\beta^{2}`` exactly.

**Composition rules that must be preserved.** The numerical
pipeline computing ``Z(\beta)`` must be:

- **algebraic**: ``+, -, \times, /, \text{power}`` all propagate
  dual numbers;
- **analytic**: `log`, `exp`, `cosh`, `sinh`, `tanh`, `sqrt`, etc.
  have dual-number overloads;
- **no comparisons on the value**: `if β > 0 then …` with a dual
  ``\beta`` compares the real part only, which is fine for
  continuous regions of ``\beta`` but breaks at discontinuities
  (none exist for our transfer-matrix path, but e.g. `abs(β)` at
  ``\beta = 0`` is non-differentiable);
- **no LAPACK routines**: `eigvals`, `eigen`, `svd`, `lu`,
  `qr`, etc. call into Fortran code that assumes `Float64` arrays,
  and the dual-number element type triggers a `MethodError`.
  The standard workaround is to avoid these entirely.

### Step 7 — Why ``\mathrm{Tr}(T^{L_{x}})``, not ``\sum\lambda_{i}^{L_{x}}``

The Ising partition function admits two equivalent forms,

```math
Z \;=\; \mathrm{Tr}(T^{L_{x}})
 \;=\; \sum_{i = 1}^{2^{L_{y}}}\lambda_{i}^{L_{x}},
```

where ``\lambda_{i}`` are the eigenvalues of the $2^{L_{y}}\times
2^{L_{y}}`` symmetric transfer matrix ``T$. The second form is more
efficient at large ``L_{x}`` (the spectrum is computed once and the
power is a scalar operation), but the eigenvalue decomposition
uses LAPACK (`DSYEV`), which **does not accept dual-number
inputs**.

The first form ``Z = \mathrm{Tr}(T^{L_{x}})`` uses only matrix
multiplication and trace, both of which are generic Julia code
with full dual-number support. Computing ``T^{L_{x}}`` via
**repeated squaring** takes ``O(\log L_{x})`` matrix
multiplications, each of which is ``O((2^{L_{y}})^{3})``. For the
typical ``L_{y} = 4, L_{x} = 4`` lattice this is ``\sim`` 10 matrix
multiplications of a ``16\times 16`` matrix — fast.

The downside is that the condition number of ``T^{L_{x}}`` grows as
``\kappa(T)^{L_{x}}`` at low temperatures (see
[`transfer-matrix-symmetric-split`](transfer-matrix-symmetric-split.md)
Step 6), so the matmul-based approach can lose precision for very
low ``T`` or very large ``L_{x}``. For the ranges used by QAtlas
(`β ≲ 5, L_{x} ≤ 6`), this is not a practical concern.

The QAtlas implementation in `src/models/classical/IsingSquare/IsingSquare.jl`
wraps the matmul path behind `fetch(IsingSquare(; J, Lx, Ly),
PartitionFunction(); β)` and exposes the AD-compatible derivatives
via `FreeEnergy`, `ThermalEntropy`, `SpecificHeat` quantities.

### Step 8 — Verification: AD vs brute-force ensemble

The test `test/verification/test_ising_ad_thermodynamics.jl`
cross-checks the AD-computed thermodynamic quantities against a
**brute-force ensemble average** over all ``2^{N}`` configurations:

```math
\langle A\rangle_{\rm brute}
 \;=\; \frac{\sum_{\{s\}} A(s)\,e^{-\beta H(s)}}
            {\sum_{\{s\}} e^{-\beta H(s)}}.
```

For the ``2\times 2``, ``2\times 3``, and ``3\times 3`` square lattices
with various ``\beta \in [0.1, 2.0]``, AD and brute force agree to

```math
|A_{\rm AD} - A_{\rm brute}|\,/\,|A_{\rm brute}| \;\le\; 10^{-14}
```

— machine precision. This confirms that:

1. The transfer-matrix construction (1)–(3) of
   [`transfer-matrix-symmetric-split`](transfer-matrix-symmetric-split.md)
   is correct.
2. The dual-number chain-rule propagation through `*`, `tr`, and
   `log` is mathematically faithful.
3. No non-differentiable operation (`abs`, `if`, LAPACK call)
   interrupts the AD pipeline.

### Step 9 — Limiting-case checks

**(i) ``\beta \to 0`` (infinite temperature).** ``Z \to 2^{N}`` and
``\ln Z \to N\ln 2``, so $\langle E\rangle = -\partial_{\beta}\ln Z
\to 0`` and ``C_{v} \to 0``. Entropy ``S \to N\ln 2$ (maximal entropy).
The AD derivatives at ``\beta = 0`` reproduce these trivial values
exactly.

**(ii) ``\beta \to \infty`` (zero temperature).** Only the ground
state contributes, ``Z \sim e^{-\beta E_{\rm gs}}`` with degeneracy
``g``; so ``\ln Z \to -\beta E_{\rm gs} + \ln g``. Then $\langle
E\rangle \to E_{\rm gs}``, ``C_{v} \to 0``, ``S \to \ln g$ (third law).
The IsingSquare AD path handles this correctly, bounded only by
the condition-number limitation noted in Step 7.

**(iii) Phase transition near ``T_{c}``.** In the thermodynamic
limit, ``C_{v}`` has a logarithmic divergence at ``T_{c}``
(see [`yang-magnetization-toeplitz`](yang-magnetization-toeplitz.md)
for ``\alpha = 0`` with a log correction; independently cross-checked
by [`ising-scaling-relations`](ising-scaling-relations.md)). On a
finite ``L_{x}\times L_{y}`` lattice the divergence is rounded into
a finite peak; the AD path reproduces the peak value and its
location consistently with ED.

**(iv) Internal consistency: fluctuation identity.** Compute
``\langle E^{2}\rangle`` and ``\langle E\rangle^{2}`` separately via
brute force, and compare to the AD-computed
``\partial^{2}\ln Z/\partial\beta^{2}``. Agreement to ``10^{-14}`` is
another cross-check of equation (2).

---

## References

- K. Huang, *Statistical Mechanics*, 2nd ed. (Wiley, 1987), Ch.
  14. Canonical ensemble + partition-function derivatives
  (equations (1)–(4)).
- J. Revels, M. Lubin, T. Papamarkou, *Forward-mode automatic
  differentiation in Julia*, arXiv:1607.07892 (2016). ForwardDiff.jl
  implementation and dual-number arithmetic.
- M. Innes, *Don't unroll adjoint: Differentiating SSA-form
  programs*, arXiv:1810.07951 (2018). Reverse-mode AD for
  thermodynamic gradients (Zygote.jl); an alternative to ForwardDiff
  for very high derivative orders.
- R. E. Wengert, *A simple automatic derivative evaluation
  program*, Commun. ACM **7**, 463 (1964). Historical origin of
  forward-mode AD via dual numbers.
- S. Huber *et al.*, *Differentiable programming tensor networks*,
  Phys. Rev. X **9**, 031041 (2019). Modern application of AD to
  transfer-matrix / tensor-network partition functions.

## Used by

- [Automatic-differentiation method](../methods/automatic-differentiation/index.md) —
  this note is the worked example of how ``\partial_{\beta}\ln Z``
  AD extracts thermodynamic quantities.
- [IsingSquare model page](../models/classical/ising-square.md) —
  `fetch(IsingSquare(; J, Lx, Ly), FreeEnergy(); β)` and its
  sibling `SpecificHeat`, `ThermalEntropy`, … dispatch through
  the transfer-matrix + AD path described in Steps 6–7.
- [`transfer-matrix-symmetric-split`](transfer-matrix-symmetric-split.md) —
  provides the symmetric ``T`` used in ``Z = \mathrm{Tr}(T^{L_{x}})``;
  AD compatibility is one of the three motivations for the
  symmetric form.
- [Ising universality class](../universalities/ising.md) — the
  IsingSquare thermodynamic values computed by this pipeline
  cross-verify the theoretical exponents (α = 0 log, …) derived
  in [`ising-scaling-relations`](ising-scaling-relations.md).
