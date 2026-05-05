# Sudden Quench of the Local Transverse Magnetisation in the TFIM

This page documents the closed-form derivation of $\langle\sigma^x_i\rangle(t)$
exposed by `MagnetizationXLocal{:quench}` for the transverse-field Ising model
(TFIM).  The quench setup follows Calabrese–Essler–Fagotti
(J. Stat. Mech. **P07016** (2012)); the underlying Bogoliubov machinery dates
to Barouch–McCoy–Dresden (Phys. Rev. A **2** (1970) 1075).

## Setup

The Hamiltonian

$$
H(h) \;=\; -J\sum_i \sigma^z_i \sigma^z_{i+1} \;-\; h \sum_i \sigma^x_i
$$

is parameterised by the transverse field $h$ at fixed $J > 0$.  A **sudden
quench** prepares the system in the ground state $|\psi_0\rangle$ of $H(h_0)$
and lets it evolve under $H(h_f)$ for $t > 0$.  The observable of interest is

$$
\boxed{\;
\langle\sigma^x_i\rangle(t)
\;=\; \langle \psi_0 | e^{i H(h_f) t}\, \sigma^x_i\, e^{-i H(h_f) t} | \psi_0 \rangle.
\;}
$$

## Open-boundary route — Majorana covariance evolution

After the Jordan–Wigner mapping discussed in
[`jw-tfim-bdg`](jw-tfim-bdg.md), $H(h)$ is quadratic in Majorana operators

$$
\gamma_{2i-1} = c_i + c_i^\dagger,\qquad \gamma_{2i} = i (c_i^\dagger - c_i),
\qquad H(h) \;=\; \tfrac{i}{4}\,\sum_{ab}\, [h]_{ab}\, \gamma_a \gamma_b,
$$

where the $2N \times 2N$ real antisymmetric matrix $[h]$ is built by
`_majorana_ham(N, J, h)`.  Its only non-zero entries are
$[h]_{2i-1,2i} = 2h$ and $[h]_{2i,2i+1} = 2J$.  In this representation
$\sigma^x_i = -i\,\gamma_{2i-1}\gamma_{2i}$.

Define the ground-state Majorana covariance

$$
\Sigma_{ab}^{(0)} \;\equiv\; -i\,\bigl(\langle \gamma_a \gamma_b\rangle_{\rm GS} - \delta_{ab}\bigr)
\;=\; -i\,\operatorname{sign}(i\,[h_0]),
$$

provided by `_majorana_covariance_gs(_majorana_ham(N, J, h_0))`.  Since
$[h_f]$ is real antisymmetric, time evolution is the orthogonal rotation

$$
\gamma_a(t) \;=\; R(t)_{ab}\,\gamma_b,
\qquad R(t) \;=\; \exp([h_f]\,t) \in \mathrm{SO}(2N).
$$

Combining these gives the time-evolved covariance

$$
\boxed{\;
\Sigma(t) \;=\; R(t)\,\Sigma^{(0)}\,R(t)^{\!T},
\;}
$$

and reading off the local magnetisation:

$$
\boxed{\;
\langle\sigma^x_i\rangle(t)
\;=\; -i\,\langle \gamma_{2i-1}(t)\,\gamma_{2i}(t)\rangle
\;=\; \Sigma(t)_{2i-1,\,2i}.
\;}
$$

This is what `_tfim_sigma_x_quench_obc` returns: one $2N \times 2N$
eigendecomposition for $\Sigma^{(0)}$ plus one matrix exponential for $R(t)$
per call.  The answer is exact at every finite $N$.

## Infinite-volume route — closed-form $k$-integral

In the thermodynamic limit translation invariance decouples the $(k, -k)$
mode pairs and the per-mode Bogoliubov rotation can be performed analytically.
With

$$
\varepsilon_k(h) \;=\; h - J \cos k,
\qquad
\Delta_k \;=\; J \sin k,
\qquad
\Lambda_k(h) \;=\; 2\sqrt{\varepsilon_k(h)^2 + \Delta_k^2},
$$

the Bogoliubov angle $\theta_k(h)$ is fixed (modulo a quadrant choice we
resolve via `atan2`) by

$$
2\,\theta_k(h) \;=\; \operatorname{atan2}\bigl(2\Delta_k,\ 2\varepsilon_k(h)\bigr),
$$

so that $\cos(2\theta_k) = 2\varepsilon_k/\Lambda_k$ and
$\sin(2\theta_k) = 2\Delta_k/\Lambda_k$.  In equilibrium at $T = 0$ the
ground-state transverse magnetisation is

$$
\langle\sigma^x\rangle_{\rm GS}(h)
\;=\; \frac{1}{\pi}\int_0^\pi \cos\bigl(2\theta_k(h)\bigr)\,dk,
$$

which agrees with the formula already coded as `MagnetizationX, Infinite`.

The post-quench state has, in the post-quench Bogoliubov basis, an excitation
density per mode pair $\sin^2(\Delta\theta_k)$ where

$$
\Delta\theta_k \;\equiv\; \theta_k(h_f) - \theta_k(h_0),
$$

and the time-dependent expectation of $\sigma^x$ within that pair is the
standard two-level Rabi-like oscillation.  Summing over $k$ gives

$$
\boxed{\;
\langle\sigma^x\rangle(t)
\;=\; \frac{1}{\pi}\int_0^\pi dk\;\Bigl[\,
\cos\bigl(2\theta_k^f\bigr)\,\cos\bigl(2\Delta\theta_k\bigr)
\;+\; \sin\bigl(2\theta_k^f\bigr)\,\sin\bigl(2\Delta\theta_k\bigr)\,
\cos\bigl(2\,\Lambda_k^f\,t\bigr)
\,\Bigr],
\;}
$$

with $\theta_k^f \equiv \theta_k(h_f)$ and $\Lambda_k^f \equiv \Lambda_k(h_f)$.
This is `_tfim_sigma_x_quench_infinite`, evaluated by adaptive Gauss–Kronrod
quadrature.

### Sanity checks

| Limit | Result |
|---|---|
| $t = 0$ | $\frac{1}{\pi}\int_0^\pi \cos(2\theta_k^f - 2\Delta\theta_k)\,dk = \frac{1}{\pi}\int_0^\pi \cos(2\theta_k(h_0))\,dk = \langle\sigma^x\rangle_{\rm GS}(h_0)$ |
| $h_0 = h_f$ | $\Delta\theta_k = 0 \Rightarrow$ time-independent $= \langle\sigma^x\rangle_{\rm GS}(h_f)$ |
| $t \to \infty$ time average | $\langle\sigma^x\rangle_{\rm GGE} = \frac{1}{\pi}\int_0^\pi \cos(2\theta_k^f)\,\cos(2\Delta\theta_k)\,dk$ (diagonal ensemble) |

The OBC route converges to the infinite-volume integral while the
**light cone** from each open boundary has not yet reached the
observation site $i$.  The relevant condition at site $i$ is

$$
t \, v_{\max} \,<\, d(i),\qquad d(i) \,=\, \min\bigl(i, \, N - i + 1\bigr),
$$

where $d(i)$ is the distance to the nearer boundary and
$v_{\max} \le 2\,\max(J, h)$ is the maximal group velocity in either
phase of the TFIM.  For the central site $i = N/2$ the bound becomes
$t\,v_{\max} < N/2$, which is the convention quoted in the
verification tests.  Outside this window the OBC value contains
boundary-reflection contributions not present in the infinite-volume
integral and a finite (controllable) deviation appears.

## API

```julia
struct MagnetizationXLocal{M} <: AbstractQuantity end          # M ∈ {:equilibrium, :quench}
MagnetizationXLocal()              = MagnetizationXLocal{:equilibrium}()
MagnetizationXLocal(:quench)       = MagnetizationXLocal{:quench}()

# OBC, single (i, t) point
fetch(model_f::TFIM, ::MagnetizationXLocal{:quench}, ::OBC;
      initial::TFIM, i::Int, t::Real, kwargs...) -> Float64

# Infinite, single t point
fetch(model_f::TFIM, ::MagnetizationXLocal{:quench}, ::Infinite;
      initial::TFIM, t::Real, kwargs...) -> Float64
```

The pre-quench TFIM `initial` and the post-quench `model_f` must share `J`;
a $J \to J'$ jump is rejected with `ArgumentError`.  Sweep callers should
hoist the Majorana eigendecomposition / quadrature setup themselves; the
public `fetch` recomputes per call.

## References

- E. Barouch, B. McCoy, M. Dresden, *Statistical Mechanics of the* XY
  *Model. I*, Phys. Rev. A **2**, 1075 (1970).
- P. Calabrese, F. H. L. Essler, M. Fagotti, *Quantum quench in the
  transverse field Ising chain: I. Time evolution of order parameter
  correlators*, J. Stat. Mech. **P07016** (2012).
- I. Peschel, *Calculation of reduced density matrices from correlation
  functions*, J. Phys. A **36**, L205 (2003).
