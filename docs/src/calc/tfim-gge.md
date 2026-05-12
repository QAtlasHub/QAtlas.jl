# Generalised Gibbs Ensemble for the TFIM Quench

## Why GGE, not Gibbs

For a generic non-integrable Hamiltonian, the eigenstate thermalisation
hypothesis predicts that *every* local observable relaxes to the
microcanonical (equivalently, canonical) Gibbs ensemble at the
temperature fixed by the initial-state energy density.  Free-fermion
models like the TFIM are **integrable** — they have an extensive
family of mutually commuting conserved charges $\{n_k\}$
(quasiparticle occupations).  The long-time average of a local
observable in a quench with these conserved charges is therefore
described by the **Generalised Gibbs Ensemble**

$$
\rho_{\rm GGE} \;=\; \frac{1}{Z}\,\exp\!\Bigl(-\sum_k \lambda_k\, n_k\Bigr),
$$

with one Lagrange multiplier $\lambda_k$ per conserved mode (Rigol et
al., PRL **98**, 050405 (2007)).  Specialising to the post-quench
TFIM, $\rho_{\rm GGE}$ becomes diagonal in the post-quench
quasiparticle basis with occupations frozen at their initial-state
values — the *diagonal ensemble*.

A single canonical Gibbs ensemble cannot match every $\langle
n_k\rangle$ simultaneously (only one parameter $\beta$), so it
generally fails to predict the long-time observables of an integrable
quench.  See Calabrese, Essler, Fagotti J. Stat. Mech. (2012) P07016 /
P07022 for a detailed TFIM analysis with non-trivial counterexamples
to canonical thermalisation.

## Setup: TFIM $h$-quench

Pre-quench Hamiltonian:

$$
H_0 \;=\; -J \sum_i \sigma^z_i\sigma^z_{i+1} \;-\; h_0 \sum_i \sigma^x_i,
$$

initial state $|\psi_0\rangle = |\text{GS}(H_0)\rangle$.  Post-quench
Hamiltonian $H_f$ is identical with $h_0 \to h_f$ (Ising coupling $J$
is held fixed).  At each momentum $k \in [0,\pi]$ the Bogoliubov
diagonalisation defines the angle

$$
2\theta_k(h) \;=\; \operatorname{atan2}\bigl(J\sin k,\; h - J\cos k\bigr),
\qquad
\Lambda_k(h) \;=\; 2\sqrt{J^2 + h^2 - 2 J h \cos k}.
$$

Because $H_0$ and $H_f$ share the momentum decomposition, every mode
occupation in the post-quench basis is conserved:

$$
\boxed{\;
n_k \;=\; \sin^2\!\bigl(\theta_k(h_0) - \theta_k(h_f)\bigr).
\;}
$$

## Closed-form GGE expectations

The diagonal-ensemble expectation of any quadratic observable in the
post-quench basis depends only on $\{n_k\}$.  For the per-site energy
and the bulk transverse magnetisation:

$$
\boxed{\;
\varepsilon_{\rm GGE} \;=\; -\frac{1}{\pi}\int_0^\pi\!\!dk\,
   \frac{\Lambda_k(h_f)}{2}\,\bigl(1 - 2 n_k\bigr),
\;}
$$

$$
\boxed{\;
\langle\sigma^x\rangle_{\rm GGE} \;=\;
   \frac{2}{\pi}\int_0^\pi\!\!dk\,
   \frac{h_f - J\cos k}{\Lambda_k(h_f)}\,
   \bigl(1 - 2 n_k\bigr).
\;}
$$

The factor $(1 - 2 n_k)$ replaces $\tanh(\beta\Lambda/2)$ of the
equilibrium expressions.  In the no-quench limit $h_0 = h_f$ one has
$n_k \equiv 0$ and the right-hand sides reduce to the $T = 0$
ground-state expressions implemented in
`TFIM.jl` / `TFIM_thermal.jl`.

## Energy conservation

The post-quench energy is a constant of motion:

$$
\langle\psi_0 \mid H_f \mid \psi_0\rangle \;=\;
\langle H_f\rangle_{\rm GGE},
$$

— the GGE energy *is* the (time-independent) initial-state expectation
of $H_f$.  This is a non-trivial cross-check: re-deriving
$\langle\psi_0|H_f|\psi_0\rangle/N$ from the raw BdG matrix elements
gives

$$
\frac{\langle\psi_0|H_f|\psi_0\rangle}{N}
   \;=\; -\frac{1}{\pi}\int_0^\pi\!\!dk\,
   \Bigl[(h_f - J\cos k)\cos(2\theta_0(k))
        + J\sin k\,\sin(2\theta_0(k))\Bigr],
$$

which is algebraically identical to the GGE form above but uses an
entirely different trigonometric branch — making it an excellent
regression test for any sign / branch error in the implementation.
The standalone test file `test/standalone/test_tfim_gge.jl` checks
this explicitly.

## API

```julia
m_0 = TFIM(J = 1.0, h = 2.0)   # pre-quench
m_f = TFIM(J = 1.0, h = 0.5)   # post-quench

# Per-site energy density of the relaxed state:
fetch(m_f, GGEValue(Energy()), Infinite(); initial = m_0)

# Stationary transverse magnetisation:
fetch(m_f, GGEValue(MagnetizationX()), Infinite(); initial = m_0)
```

The `initial::TFIM` keyword is required.  Mismatched Ising couplings
(`m_0.J != m_f.J`) raise a `DomainError` — only $h$-quenches are
covered by the closed forms above.

## References

- M. Rigol, V. Dunjko, V. Yurovsky, M. Olshanii, *Relaxation in a
  Completely Integrable Many-Body Quantum System: An Ab Initio Study
  of the Dynamics of the Highly Excited States of 1D Lattice Hard-Core
  Bosons*, **PRL 98**, 050405 (2007).
- P. Calabrese, F. H. L. Essler, M. Fagotti, *Quantum Quench in the
  Transverse Field Ising Chain*, J. Stat. Mech. (2012) **P07016**,
  **P07022**.
- M. Fagotti, F. H. L. Essler, *Reduced density matrix after a quantum
  quench*, **PRB 87**, 245107 (2013).
