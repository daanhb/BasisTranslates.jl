
"""
A dictionary of periodic B-splines on the interval [0,1].
"""
struct PeriodicBSplines{T} <: BasisTranslates.PeriodicTranslates{T,T}
    degree   ::  Int
    n       ::  Int
end

PeriodicBSplines(; degree, n) = PeriodicBSplines(n, degree)

PeriodicBSplines(degree::Int, n::Int) =
    PeriodicBSplines{Float64}(degree, n)

BasisFunctions.size(Φ::PeriodicBSplines) = (Φ.n,)

spline_degree(Φ::PeriodicBSplines) = Φ.degree
spline_order(Φ) = spline_degree(Φ)+1


# We use as a kernel the BSpline on [0,degree+1]
BasisTranslates.kernel_support(Φ::PeriodicBSplines{T}) where {T} =
    Interval(zero(T), spline_degree(Φ)+one(T))

# The basis functions are centered (and dilated) splines.
# Exception for even degree splines: we align the origin with a node,
# rather than the origin being in between two nodes. This allows for
# exact conversion between spline spaces of different degree.
semicentered_shift(degree::Int, ::Type{T}) where {T} =
    iseven(degree) ? T(degree)/2 : T(degree+1)/2
semicentered_shift(Φ::PeriodicBSplines{T}) where {T} =
    semicentered_shift(spline_degree(Φ), T)

# the map from [0,1] to the kernel domain takes into account the
# shift by i, the shift of the kernel, and a scaling by N
BasisTranslates.translate_map(Φ::PeriodicBSplines, i) =
    AffineMap(length(Φ), -(i-1) + semicentered_shift(Φ))

# this version is slightly more accurate than the default
function BasisFunctions.support(Φ::PeriodicBSplines, idx)
    shift = semicentered_shift(Φ)
    subdomain = Interval(-shift+(idx-1), -shift+(idx-1)+spline_degree(Φ)+1)/length(Φ)
    PeriodicInterval(subdomain, support(Φ))
end

BasisTranslates.kernel_eval(Φ::PeriodicBSplines{T}, x) where {T} =
    eval_bspline(spline_degree(Φ), x, T)

BasisTranslates.kernel_eval_derivative(Φ::PeriodicBSplines{T}, x, order) where {T} =
    eval_bspline_derivative(spline_degree(Φ), x, order, T)
