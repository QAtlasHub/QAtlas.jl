# ─────────────────────────────────────────────────────────────────────────────
# XXZ1D — Exact two-spinon longitudinal dynamical structure factor (Pérez Castillo 2020)
#
# Implements the exact two-spinon longitudinal dynamical structure factor
# S^{zz}(Q, ω) of the spin-1/2 antiferromagnetic XXZ chain in the massive
# regime Δ > 1 (thermodynamic limit, Infinite boundary conditions),
# as derived in Isaac Pérez Castillo, arXiv:2005.10729.
# ─────────────────────────────────────────────────────────────────────────────

using QuadGK: quadgk

function _xxz_elliptic_params(Δ::Real)
    ε = acosh(float(Δ))
    if ε > 1.0
        # Direct summation for large ε
        qe = exp(-ε)

        theta2 = 0.0
        for n in 0:12
            theta2 += qe^(n * (n + 1))
        end
        theta2 *= 2.0 * qe^(0.25)

        theta3 = 1.0
        for n in 1:12
            theta3 += 2.0 * qe^(n^2)
        end

        theta4 = 1.0
        for n in 1:12
            theta4 += 2.0 * (-1.0)^n * qe^(n^2)
        end

        k = (theta2 / theta3)^2
        kp = (theta4 / theta3)^2
        K_val = (pi / 2) * theta3^2
        Kp_val = (ε / pi) * K_val
        return qe, k, kp, K_val, Kp_val, ε
    else
        # Modular transformation for small ε (Δ -> 1+) to ensure numerical stability
        qp = exp(-pi^2 / ε)

        theta2 = 0.0
        for n in 0:8
            theta2 += qp^(n * (n + 1))
        end
        theta2 *= 2.0 * qp^(0.25)

        theta3 = 1.0
        for n in 1:8
            theta3 += 2.0 * qp^(n^2)
        end

        theta4 = 1.0
        for n in 1:8
            theta4 += 2.0 * (-1.0)^n * qp^(n^2)
        end

        k = (theta4 / theta3)^2
        kp = (theta2 / theta3)^2
        Kp_val = (pi / 2) * theta3^2
        K_val = (pi / ε) * Kp_val
        qe = exp(-ε)
        return qe, k, kp, K_val, Kp_val, ε
    end
end

function _incomplete_elliptic_F(ϕ::Real, k::Real)
    integrand(θ) = 1.0 / sqrt(1.0 - k^2 * sin(θ)^2)
    val, _ = quadgk(integrand, 0.0, ϕ; atol=1e-14, rtol=1e-14)
    return val
end

function _neville_theta_d(u::Real, K_val::Real, qe::Real)
    num = 1.0
    for n in 1:12
        num += 2.0 * qe^(n^2) * cos(n * pi * u / K_val)
    end

    den = 1.0
    for n in 1:12
        den += 2.0 * qe^(n^2)
    end

    return num / den
end

function _thetaA_sq(β::Real, K_val::Real, ε::Real)
    sum_val = 0.0
    # Dynamic number of terms to ensure convergence as ε -> 0
    max_k = max(100, round(Int, 10.0 / ε))
    for k in 1:max_k
        theta = 2.0 * pi * k * β / K_val
        x = k * ε
        # Numerically stable form avoiding exp(x) overflow and catastrophic cancellation
        num = cos(theta) * expm1(-2x)^2 - 4.0 * exp(-2x) * sin(theta / 2)^2
        den = -expm1(-4x) * (1.0 + exp(-2x))

        ratio = 2.0 * num / den
        u_k = ratio / k
        sub_k = 2.0 * cos(theta) / k
        sum_val += (u_k - sub_k)
    end

    val_exp = exp(-sum_val)
    sin_term = 4.0 * sin(pi * β / K_val)^2
    return sin_term * val_exp
end

function _xxz_szz_term(
    Q_val::Real,
    ω::Real,
    qe::Real,
    k::Real,
    kp::Real,
    K_val::Real,
    ε::Real,
    κ::Real,
    Δ::Real,
    σ::Real,
)
    sinQ = sin(Q_val)
    if abs(sinQ) < 1e-12
        return 0.0
    end

    # In units where J=1
    I_val = (K_val / pi) * sqrt(Δ^2 - 1.0)
    ω0 = (2.0 * I_val / (1.0 + κ)) * sinQ

    B = sqrt(ω^2 - κ^2 * ω0^2) * sqrt(ω^2 - ω0^2)

    W_arg = κ^2 * (ω0^4 / ω^4) - (B / ω^2 + cos(Q_val))^2
    if W_arg <= 0.0
        return 0.0
    end
    W = sqrt(W_arg)

    dn_arg =
        ((1.0 + cos(Q_val)) / abs(sinQ)) *
        sqrt((ω^2 - κ^2 * ω0^2 + B) / (ω^2 + κ^2 * ω0^2 - B))
    dn_arg_clipped = clamp(dn_arg, kp, 1.0)

    phi_arg = sqrt(1.0 - dn_arg_clipped^2) / k
    phi = asin(clamp(phi_arg, 0.0, 1.0))

    beta = _incomplete_elliptic_F(phi, k)

    theta_d_val = _neville_theta_d(beta, K_val, qe)
    thetaA_val = _thetaA_sq(beta, K_val, ε)

    num = (1.0 + cos(Q_val)) * thetaA_val * (ω^2 - σ * (B - κ * ω0^2))
    den = W * (theta_d_val^2) * (Δ - σ * cos(pi * beta / K_val))

    return num / den
end

"""
    fetch(model::XXZ1D, ::ZZStructureFactor, ::Infinite;
          q::Real, ω::Real, J::Real = model.J, method::Symbol = :exact_2spinon,
          kwargs...) -> Float64

Exact two-spinon longitudinal dynamical structure factor S^{zz}(q, ω) of the
spin-1/2 antiferromagnetic XXZ chain in the massive regime Δ > 1.
"""
function fetch(
    model::XXZ1D,
    ::ZZStructureFactor,
    ::Infinite;
    q::Real,
    ω::Real,
    J::Real=model.J,
    method::Symbol=:exact_2spinon,
    kwargs...,
)
    J > 0 || throw(DomainError(J, "XXZ1D ZZStructureFactor requires J > 0; got J = $J."))
    model.Δ > 1.0 || throw(
        DomainError(
            model.Δ,
            "XXZ1D exact 2-spinon longitudinal DSF is only defined in the massive regime Δ > 1.",
        ),
    )

    if method !== :exact_2spinon
        throw(
            ArgumentError(
                "XXZ1D ZZStructureFactor Infinite: unknown method :$method; supported is :exact_2spinon.",
            ),
        )
    end

    # Scale ω by J to compute in units J=1
    ωf = float(ω)
    Jf = float(J)
    ω_unit = ωf / Jf

    # Parameters in units J=1
    qe, k, kp, K_val, Kp_val, ε = _xxz_elliptic_params(model.Δ)
    κ = (1.0 - kp) / (1.0 + kp)
    Q_κ = acos(κ)

    Q_norm = mod(float(q), 2*pi)
    if Q_norm > pi
        Q_norm = 2*pi - Q_norm
    end

    I_val = (K_val / pi) * sqrt(model.Δ^2 - 1.0)

    omega0(q_val) = (2.0 * I_val / (1.0 + κ)) * sin(q_val)
    omegaplus(q_val) = (2.0 * I_val / (1.0 + κ)) * sqrt(1.0 + κ^2 + 2.0 * κ * cos(q_val))
    omegaminus(q_val) = (2.0 * I_val / (1.0 + κ)) * sqrt(1.0 + κ^2 - 2.0 * κ * cos(q_val))

    in_C_plus = false
    if Q_norm >= Q_κ
        lo = (Q_norm <= pi - Q_κ) ? omega0(Q_norm) : omegaplus(Q_norm)
        up = omegaminus(Q_norm)
        if lo <= ω_unit <= up
            in_C_plus = true
        end
    end

    in_C_minus = false
    Q_reflected = pi - Q_norm
    if Q_reflected >= Q_κ
        lo = (Q_reflected <= pi - Q_κ) ? omega0(Q_reflected) : omegaplus(Q_reflected)
        up = omegaminus(Q_reflected)
        if lo <= ω_unit <= up
            in_C_minus = true
        end
    end

    if !in_C_plus && !in_C_minus
        return 0.0
    end

    B = sqrt(ω_unit^2 - κ^2 * omega0(Q_norm)^2) * sqrt(ω_unit^2 - omega0(Q_norm)^2)
    prefactor = sqrt(qe) * k * (ω_unit^2 + κ * omega0(Q_norm)^2 + B) / (ω_unit^3 * B)

    term_plus = 0.0
    if in_C_plus
        term_plus =
            ((1.0 - κ) / (1.0 + κ)) *
            _xxz_szz_term(Q_norm, ω_unit, qe, k, kp, K_val, ε, κ, model.Δ, 1.0)
    end

    term_minus = 0.0
    if in_C_minus
        term_minus = _xxz_szz_term(
            Q_reflected, ω_unit, qe, k, kp, K_val, ε, κ, model.Δ, -1.0
        )
    end

    S_unit = prefactor * (term_plus + term_minus)
    # Divide by 4.0 to convert from Pauli-spin operator convention to standard Spin S=1/2 representation
    return (S_unit / Jf) / 4.0
end
