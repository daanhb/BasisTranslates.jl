
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


"""
Refinable function as defined in the 1992 paper by Cohen, Daubechies and
Feauveau "Biorthogonal bases of compactly supported wavelets".

In this package the convention is that the primal function is the B-Spline of
order `p`. The dual function of order `q` is just one of the several
possibilities.
"""
struct CDFDual{T} <: Refinable{T}
    p               ::  Int
    q               ::  Int
    coefficients    ::  VectorSequence{T}
end

refinable_coeff(φ::CDFDual) = φ.coefficients

function cdfdual_filter(p::Int, q::Int, ::Type{T} = Float64) where T
    z = [-1]
    if p == 1
        if q == 1
            z = [1, 1]
        elseif q == 3
            z = [-1, 1, 8, 8, 1, -1]
        elseif q == 5
            z = [3, -3, -22, 22, 128, 128, 22, -22, -3, 3]
        end
    elseif p == 2
        if q == 2
            z = [-1, 2, 6, 2, -1]
        elseif q == 4
            z = [3, -6, -16, 38, 90, 38, -16, -6, 3]
        elseif q == 6
            z = [-5, 10, 34, -78, -123, 324, 700, 324, -123, -78, 34, 10, -5]
        end
    elseif p == 3
        if q == 1
            z = [-1, 3, 3, -1]
        elseif q == 3
            z = [3, -9, -7, 45, 45, -7, -9, 3]
        elseif q == 5
            z = [-5, 15, 19, -97, -26, 350, 350, -26, -97, 19, 15, -5]
        end
    elseif p == 4
        if q == 2
            z = [3, -12, 5, 40, 5, -12, 3]
        elseif q == 4
            z = [-10, 40, -2, -192, 140, 560, 140, -192, -2, 40, -10]
        elseif q == 6
            z = [35, -140, -55, 920, -557, -2932, 2625, 8400, 2625, -2932, -557, 920, -55, -140, 35]
        end
    elseif p == 5
        if q == 1
            z = [3, -15, 20, 20, -15, 3]
        elseif q == 3
            z = [-5, 25, -26, -70, 140, 140, -70, -26, 25, -5]
        elseif q == 5
            z = [35, -175, 120, 800, -1357, -1575, 4200, 4200, -1575, -1357, 800, 120, -175, 35]
        end
    elseif p == 6
        if q == 2
            z = [-5, 30, -56, -14, 154, -14, -56, 30, -5]
        elseif q == 4
            z = [35, -210, 330, 470, -1827, 252, 3948, 252, -1827, 470, 330, -210, 35]
        elseif q == 6
            z = [-63, 378, -476, -1554, 4404, 1114, -13860, 4158, 28182, 4158, -13860, 1114, 4404, -1554, -476, 378, -63]
        end
    end
    if z == [-1]
        error("Combination (p,q) = ($(p), $(q)) not implemented for CDF duals.")
    end
    # make the coefficients of type T and ensure they sum to 2
    data = (z * sqrt(T(2))) / sum(z)
    CompactSequence(data)
end

function cdfdual(p::Int, q::Int, ::Type{T} = Float64) where T
    coef = cdfdual_filter(p, q, T)
    CDFDual(p, q, coef)
end
