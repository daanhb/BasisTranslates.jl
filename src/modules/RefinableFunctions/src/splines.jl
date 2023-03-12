
using BasisTranslates.BSplines: BSpline

"Refinable coefficients of a bspline."
bspline_refinable_coefficients(degree::Int, ::Type{T} = Float64) where {T} =
    [binomial(degree+1,k) for k in 0:degree+1] / 2^(degree+one(T)/2)

"Causal filter associated with a bspline of the given degree."
function bspline_causal_filter(degree::Int, ::Type{T} = Float64) where {T}
    coef = bspline_refinable_coefficients(degree, T)
    VectorSequence(coef, 0:degree+1)
end

isrefinable(φ::BSpline) = true
refinable_coeff(φ::BSpline) = bspline_causal_filter(spline_degree(φ), numtype(φ))
refinable_moment(φ::BSpline) = one(numtype(φ))



struct CDFDual{T} <: Refinable{T}
    p               ::  Int
    q               ::  Int
    coefficients    ::  VectorSequence{T}
end

refinable_coeff(φ::CDFDual) = φ.coefficients
