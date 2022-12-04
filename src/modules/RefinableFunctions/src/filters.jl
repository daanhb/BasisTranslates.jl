
function db2(::Type{T} = Float64) where {T}
    sq2 = sqrt(T(2))
    sq3 = sqrt(T(3))
    coef = 1/sq2*[(1+sq3)/4, (3+sq3)/4, (3-sq3)/4, (1-sq3)/4]
    VectorSequence(coef)
end

bspline_refinable_coefficients(degree::Int, ::Type{T} = Float64) where {T} =
    [binomial(degree+1,k) for k in 0:degree+1] / 2^(degree+one(T)/2)

function bspline(degree::Int, ::Type{T} = Float64) where {T}
    coef = bspline_refinable_coefficients(degree, T)
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
