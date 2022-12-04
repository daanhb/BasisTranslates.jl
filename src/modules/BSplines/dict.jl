
using BasisTranslates:
    RegularGrid

import BasisTranslates:
    kernel,
    parent_kernel,
    kernel_support,
    centers,
    map_to_kernel

"Supertype of BSpline translates."
abstract type BSplineTranslates{T,EXT} <: BasisTranslates.Translates{T,T,EXT} end

"The degree of the spline functions."
spline_degree(Φ::BSplineTranslates) = Φ.degree

"The order of a spline is its degree plus one."
spline_order(Φ::BSplineTranslates) = spline_degree(Φ)+1


"""
A dictionary of B-splines with equispaced nodes.
"""
struct RegularBSplines{T} <: BSplineTranslates{T,:simple}
    nodes   ::  RegularGrid{T}
    degree  ::  Int
end

RegularBSplines(nodes::AbstractVector; degree = 1) =
    RegularBSplines(nodes, degree)

"The nodes of regular B-splines."
nodes(Φ::RegularBSplines) = Φ.nodes

kernel(Φ::RegularBSplines{T}) where {T} = CenteredBSpline{T}(spline_degree(Φ))

Base.step(Φ::RegularBSplines) = step(Φ.nodes)

centers(Φ::RegularBSplines) = nodes(Φ)

# Implement a linear scaling
function map_to_kernel(Φ::RegularBSplines, i)
    h = step(Φ)
    c = center(Φ, i)
    AffineMap(inv(h), -c/h)
end


"""
A dictionary of periodic B-splines on the interval `[0,1]`.
"""
struct PeriodicBSplines{T} <: BSplineTranslates{T,:unitperiodic}
    n       ::  Int
    degree  ::  Int
end

PeriodicBSplines(n::Int; degree) = PeriodicBSplines{Float64}(n; degree)
PeriodicBSplines{T}(n::Int; degree::Int) where {T} =
    PeriodicBSplines{T}(n, degree)

Base.size(Φ::PeriodicBSplines) = (Φ.n,)

spline_degree(Φ::PeriodicBSplines) = Φ.degree

parent_kernel(Φ::PeriodicBSplines{T}) where {T} = CenteredBSpline{T}(Φ.degree)
