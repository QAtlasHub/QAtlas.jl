# Loschmidt Echo and Dynamical Quantum Phase Transitions in the TFIM

## Setting

After preparing the system in the ground state ``|\psi_0\rangle`` of an
"initial" Hamiltonian
``H_0 = -J\sum_i \sigma^z_i \sigma^z_{i+1} - h_0 \sum_i \sigma^x_i``
and **suddenly quenching** the field to ``h_f``, the time-evolved state
``|\psi(t)\rangle = e^{-i H_f t}|\psi_0\rangle`` has overlap with the
initial state given by the **Loschmidt amplitude**

```math
G(t) \;=\; \langle \psi_0 | e^{-i H_f t} | \psi_0 \rangle.
```

The **Loschmidt echo** is its squared modulus and the
**rate function** (a quench dynamical free-energy density) is

```math
L(t) \;=\; |G(t)|^2,
\qquad
\lambda(t) \;=\; -\frac{1}{N} \log L(t).
```

Heyl, Polkovnikov, and Kehrein observed [^Heyl2013] that
``\lambda(t)`` exhibits **non-analytic cusps** at a discrete set of
critical times ``t_n^*`` whenever the quench crosses the equilibrium
quantum critical point — these are **dynamical quantum phase
transitions** (DQPTs). The TFIM is the canonical solvable example;
the review [^Heyl2018] traces every subsequent generalization back to
this case.

## Per-mode product structure

Both ``H_0`` and ``H_f`` share the Jordan–Wigner / Bogoliubov
decoupling into independent fermionic modes. Mode ``n`` in ``H_0``
rotates onto mode ``n`` in ``H_f`` by a Bogoliubov angle difference
``\Delta\theta_n = \theta_n^{(0)} - \theta_n^{(f)}``, and only this
two-dimensional subspace contributes to the per-mode amplitude. The
full echo factorizes over modes,

```math
\boxed{\;
L(t) \;=\; \prod_n \Bigl|
\cos^2\!\Delta\theta_n
+ \sin^2\!\Delta\theta_n \, e^{-2 i \Lambda_n^{(f)} t}
\Bigr|^2,
\;}
```

with ``\Lambda_n^{(f)}`` the H_f single-quasiparticle energy of mode ``n``.

## Infinite (continuous-``k``) form

In the thermodynamic limit, momentum becomes continuous and the
product turns into an integral. With

```math
\Lambda_k(h) = 2\sqrt{J^2 + h^2 - 2 J h \cos k},
\qquad
\tan(2\theta_k(h)) = \frac{J \sin k}{h - J \cos k},
```

the rate function reads [^Heyl2013, eq. (3)]

```math
\boxed{\;
\lambda(t) \;=\; -\frac{1}{2\pi}\int_0^\pi
\log\!\Bigl|
\cos^2\!\Delta\theta_k
+ \sin^2\!\Delta\theta_k \, e^{-2 i \Lambda_k(h_f) t}
\Bigr|^2 \, dk.
\;}
```

QAtlas evaluates this integral with `QuadGK.quadgk` directly; the
DQPT cusp is an integrable logarithmic singularity that adaptive
Gauss–Kronrod handles without special bookkeeping.

## DQPT critical times

Set ``\alpha_k = \cos^2\Delta\theta_k``, ``\beta_k = \sin^2\Delta\theta_k``.
The integrand ``|\alpha_k + \beta_k e^{-2 i \Lambda_f t}|^2`` vanishes
iff ``\alpha_k = \beta_k`` **and** ``e^{-2 i \Lambda_f t} = -1``, i.e.

```math
\cos(2\Delta\theta_{k^*}) = 0,
\qquad
2 \Lambda_{k^*}^{(f)} t_n^* = (2n+1)\pi.
```

So the DQPT critical times are

```math
\boxed{\;
t_n^* \;=\; \frac{\pi (n + \tfrac{1}{2})}{\Lambda_{k^*}^{(f)}},
\qquad n = 0, 1, 2, \dots
\;}
```

The mode ``k^*`` exists only when ``h_0`` and ``h_f`` lie on **opposite
sides** of the QCP ``h = J`` — equivalently, the quench must traverse
the equilibrium critical point. For ``h_0, h_f`` in the same phase,
``\cos(2\Delta\theta_k)`` has no zero on ``(0, \pi)`` and ``\lambda(t)`` is
real-analytic for all ``t``.

## OBC (finite ``N``) implementation

For an open chain, the Bogoliubov diagonalization yields per-site
amplitudes ``(\phi_n, \psi_n) \in \mathbb{R}^N`` via the
Lieb–Schultz–Mattis trick

```math
(A - B)(A + B) \psi_n = \Lambda_n^2 \psi_n,
\qquad
\phi_n = (A + B) \psi_n / \Lambda_n,
```

with the bipartite blocks ``A_{ii} = 2h``, ``A_{i,i \pm 1} = -J``,
``B_{i,i+1} = +J``, ``B_{i+1,i} = -J``. The Bogoliubov rotations
``G = (\phi + \psi)/2``, ``H = (\phi - \psi)/2`` map fermion operators to
quasiparticles. The overlap of the H_0 and H_f Bogoliubov vacua
factorizes into per-mode terms whose row norms

```math
\cos^2\theta_n = \sum_m |P^{(+)}_{n,m}|^2,
\qquad
\sin^2\theta_n = \sum_m |P^{(-)}_{n,m}|^2,
```

with

```math
P^{(+)} = G^{(f)\top} G^{(0)} + H^{(f)\top} H^{(0)},
\qquad
P^{(-)} = G^{(f)\top} H^{(0)} + H^{(f)\top} G^{(0)},
```

reproduce the ``\cos^2\Delta\theta + \sin^2\Delta\theta e^{-2 i \Lambda t}``
factor in the diagonal-pair limit. For OBC, residual mode mixing is
absorbed into per-mode ``\theta_n`` consistent with unitarity
(``\cos^2\theta_n + \sin^2\theta_n = 1`` enforced by row normalization
to suppress round-off drift).

## Public API

```julia
using QAtlas

m_0 = TFIM(J = 1.0, h = 2.0)   # paramagnetic
m_f = TFIM(J = 1.0, h = 0.5)   # ferromagnetic — quench crosses h = J

# Infinite (continuous-k integral, exact in N → ∞)
λ_inf = QAtlas.fetch(m_f, LoschmidtRateFunction(), Infinite();
                     initial = m_0, t = 1.0)

# OBC finite N (BdG diagonalisation)
λ_obc = QAtlas.fetch(m_f, LoschmidtRateFunction(), OBC(64);
                     initial = m_0, t = 1.0)
L_obc = QAtlas.fetch(m_f, LoschmidtAmplitude(), OBC(64);
                     initial = m_0, t = 1.0)
```

The pre-quench Hamiltonian is passed via the `initial::TFIM` keyword;
`initial.J` must match `model_f.J` (only `h` is quenched).

## Verification

`test/standalone/test_tfim_loschmidt.jl` exercises:

1. ``L(0) = 1``, ``\lambda(0) = 0``.
2. No-quench identity: ``h_0 = h_f \implies L(t) \equiv 1, \lambda(t) \equiv 0``.
3. ``L(t) \in [0, 1]``, ``\lambda(t) \geq 0`` across a sweep.
4. DQPT cusp: ``h_0 = 2 \to h_f = 0.5`` shows a slope sign-flip and a
   local maximum at ``t_c = \pi/(2 \Lambda_{k^*}^{(f)})``, with
   left/right slope discontinuity dwarfing a smooth-region control
   point.
5. OBC large-``N`` convergence to Infinite at off-cusp ``t``.
6. Pinned reference: $\lambda(t = 1; h_0=2, h_f=0.5) \approx
   0.31693310885932685$ at `atol = 1e-8`.

## References

[^Heyl2013]: M. Heyl, A. Polkovnikov, S. Kehrein,
    *Dynamical Quantum Phase Transitions in the Transverse-Field Ising Model*,
    Phys. Rev. Lett. **110**, 135704 (2013).

[^Heyl2018]: M. Heyl,
    *Dynamical quantum phase transitions: a review*,
    Rep. Prog. Phys. **81**, 054001 (2018).
