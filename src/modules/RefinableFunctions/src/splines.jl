
using BasisTranslates.BSplines: BSpline

isrefinable(φ::BSpline) = true

"Refinable coefficients of a bspline."
bspline_refinable_refinable_coeff(degree::Int, ::Type{T} = Float64) where {T} =
    [binomial(degree+1,k) for k in 0:degree+1] / 2^(degree+one(T)/2)

"Filter associated with a bspline of the given degree."
function bspline_filter(degree::Int, ::Type{T} = Float64) where {T}
    coef = bspline_refinable_refinable_coeff(degree, T)
    if isodd(length(coef))
        # we make the basis functions symmetric
        L = length(coef)>>1
        I = -L:L
    else
        # we make the filter causal
        I = 0:length(coef)-1
    end
    VectorSequence(coef, I)
end

refinable_coeff(φ::BSpline{T}) where T = bspline_filter(φ.degree, T)
