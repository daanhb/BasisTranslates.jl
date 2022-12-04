
using BasisTranslates:
    RegularGrid

import BasisTranslates:
    kernel,
    parent_kernel,
    kernel_support,
    centers,
    map_to_kernel

"""
A dictionary of B-splines with equispaced nodes.
"""
struct RegularBSplines{T} <: BasisTranslates.Translates{T,T}
    nodes   ::  RegularGrid{T}
    degree  ::  Int
end

RegularBSplines(nodes::AbstractVector; degree = 1) =
    RegularBSplines(nodes, degree)

"The degree of the spline functions."
spline_degree(Φ::RegularBSplines) = Φ.degree

"The order of a spline is its degree plus one."
spline_order(Φ) = spline_degree(Φ)+1

"The nodes of regular B-splines."
nodes(Φ::RegularBSplines) = Φ.nodes

Base.step(Φ::RegularBSplines) = step(Φ.nodes)

centers(Φ::RegularBSplines) = nodes(Φ)

kernel(Φ::RegularBSplines{T}) where {T} = CenteredBSpline{T}(Φ.degree)

# Implement a linear scaling
function map_to_kernel(Φ::RegularBSplines, i)
    h = step(Φ)
    c = center(Φ, i)
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

parent_kernel(Φ::PeriodicBSplines{T}) where {T} = CenteredBSpline{T}(Φ.degree)
