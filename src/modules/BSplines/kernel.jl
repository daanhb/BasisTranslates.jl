
import BasisTranslates:
    kernel_eval,
    kernel_eval_derivative,
    kernel_support


"Supertype of B-spline kernels."
abstract type BSplineKernel{T} <: BasisTranslates.Kernel
end

"Kernel representing a B-spline of a certain degree."
struct BSpline{T} <: BSplineKernel{T}
    degree  ::  Int
end
BSpline(degree::Int) = BSpline{Float64}(degree)

"Kernel representing a centered B-spline of a certain degree."
struct CenteredBSpline{T} <: BSplineKernel{T}
    degree  ::  Int
end
CenteredBSpline(degree::Int) = CenteredBSpline{Float64}(degree)

function kernel_support(kern::BSpline{T}) where {T}
    a, b = support_bspline(kern.degree, T)
    a..b
end

function kernel_support(kern::CenteredBSpline{T}) where {T}
    a, b = support_centered_bspline(kern.degree, T)
    a..b
end

kernel_eval(kern::BSpline{T}, x) where {T} =
    eval_bspline(kern.degree, x, T)

kernel_eval(kern::CenteredBSpline{T}, x) where {T} =
    eval_centered_bspline(kern.degree, x, T)

kernel_eval_derivative(kern::BSpline{T}, x, order::Int) where {T} =
    eval_bspline_derivative(kern.degree, x, order, T)

kernel_eval_derivative(kern::CenteredBSpline{T}, x, order::Int) where {T} =
    eval_centered_bspline_derivative(kern.degree, x, order, T)
