
"""
A dictionary of partial fractions with simple poles.
"""
struct PartialFractions{T} <: BasisTranslates.SimpleTranslates{T,T}
    poles   ::  Vector{T}
    support ::  Domain{T}
end

PartialFractions(poles::Vector{T}) where {T} =
    PartialFractions{T}(poles, UnitInterval{T}())

poles(Φ::PartialFractions) = Φ.poles

BasisTranslates.centers(Φ::PartialFractions) = Φ.poles

BasisFunctions.support(Φ::PartialFractions) = Φ.support

kernel(Φ::PartialFractions) = PartialFractionKernel()

function BasisFunctions.dict_norm(Φ::PartialFractions, idx, p::Real)
    @boundscheck checkbounds(Φ, idx)
    @assert (p == Inf) || (p==2)
    @assert support(Φ) == 0..1
    c = -Φ.poles[idx]
    if p == Inf
        abs(1/c)
    else
        sqrt(-1/(c+1) + 1/c)
    end
end

"Evaluate a linear combination of partial fractions in the point `x`."
pf_eval(z, w, x) = sum(w[k]/(x-z[k]) for k in 1:length(z))

BasisFunctions.unsafe_eval_expansion(Φ::PartialFractions, coefficients, x) =
    pf_eval(Φ.poles, coefficients, x)
