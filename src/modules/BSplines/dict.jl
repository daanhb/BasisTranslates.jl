
using BasisTranslates:
    RegularGrid

import BasisTranslates:
    kernel,
    parent_kernel,
    kernel_support,
    translates_grid,
    translate_map

"""
A dictionary of B-splines with equispaced nodes.
"""
struct RegularBSplines{T} <: BasisTranslates.Translates{T,T}
    nodes   ::  RegularGrid{T}
    degree  ::  Int
end

RegularBSplines(nodes::AbstractVector; degree = 1) =
    RegularBSplines(nodes, degree)

spline_degree(Φ::RegularBSplines) = Φ.degree
spline_order(Φ) = spline_degree(Φ)+1
nodes(Φ::RegularBSplines) = Φ.nodes

Base.step(Φ::RegularBSplines) = step(Φ.nodes)

translates_grid(Φ::RegularBSplines) = nodes(Φ)

kernel(Φ::RegularBSplines{T}) where {T} =
    CenteredBSpline{T}(Φ.degree)

function translate_map(Φ::RegularBSplines, i)
    h = step(Φ)
    c = translate_center(Φ, i)
    AffineMap(inv(h), -c/h)
end


"""
A dictionary of periodic B-splines on the interval `[0,1]`.
"""
struct PeriodicBSplines{T} <: BasisTranslates.PeriodicTranslates{T,T}
    n       ::  Int
    degree  ::  Int
end

PeriodicBSplines(n::Int; degree) = PeriodicBSplines(n, degree)

PeriodicBSplines(n::Int, degree::Int) =
    PeriodicBSplines{Float64}(n, degree)

Base.size(Φ::PeriodicBSplines) = (Φ.n,)

spline_degree(Φ::PeriodicBSplines) = Φ.degree

parent_kernel(Φ::PeriodicBSplines{T}) where {T} =
    BSpline{T}(Φ.degree)

# We use as a kernel the BSpline on [0,degree+1]
kernel_support(Φ::PeriodicBSplines{T}) where {T} =
    Interval(zero(T), spline_order(Φ))

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
translate_map(Φ::PeriodicBSplines, i) =
    AffineMap(length(Φ), -(i-1) + semicentered_shift(Φ))

# this version is slightly more accurate than the default
function BasisFunctions.support(Φ::PeriodicBSplines, idx)
    shift = semicentered_shift(Φ)
    subdomain = Interval(-shift+(idx-1), -shift+(idx-1)+spline_order(Φ))/length(Φ)
    PeriodicInterval(subdomain, support(Φ))
end
