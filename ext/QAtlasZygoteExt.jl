# Reverse-mode AD backend for QAtlas's derived inputs (core/derivative.jl).
#
# Provided for cross-checking a ForwardDiff result with an independent
# differentiation mode — the same "two routes must agree" discipline the atlas
# applies to physical values.
#
# Expect a lower success rate than ForwardDiff on this atlas: reverse mode
# needs an adjoint for everything it traverses (dense ED, quadrature, the NLIE
# Newton solves), and for a scalar → scalar derivative it buys nothing that
# forward mode does not already give.  `preferred_backend` therefore ranks it
# BELOW ForwardDiff; a hub where Zygote throws is a Zygote coverage gap, not a
# physics finding.
module QAtlasZygoteExt

using QAtlas: QAtlas, ZygoteBackend
using Zygote: Zygote

function QAtlas.derivative(::ZygoteBackend, f, x::Real)
    g = only(Zygote.gradient(f, float(x)))
    g === nothing && error(
        "QAtlas.derivative(ZygoteBackend, …): Zygote returned no gradient at x = $x — " *
        "the fetch is not reverse-differentiable (no adjoint along its path). Use " *
        "ForwardDiffBackend() or FiniteDifference() for this hub.",
    )
    return g
end
QAtlas.backend_available(::ZygoteBackend) = true

end
