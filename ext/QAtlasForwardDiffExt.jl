# Forward-mode AD backend for QAtlas's derived inputs (core/derivative.jl).
#
# Loading ForwardDiff is all it takes to upgrade every derivative-supplied
# constraint edge from the finite-difference fallback to machine precision —
# and, via `default_rtol`, to the tighter tolerance that accuracy earns.
#
# Forward mode is the right mode here: the derived inputs are scalar → scalar
# (`dF/dT`, `dβF/dβ`, `dM/dh`), and it needs nothing from a `fetch` beyond
# being generic in its argument type, which 273 of the 304 `beta::` annotations
# in src/ already are.
module QAtlasForwardDiffExt

using QAtlas: QAtlas, ForwardDiffBackend
using ForwardDiff: ForwardDiff

QAtlas.derivative(::ForwardDiffBackend, f, x::Real) = ForwardDiff.derivative(f, float(x))
QAtlas.backend_available(::ForwardDiffBackend) = true

end
