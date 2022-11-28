
"""
A dictionary of partial fractions with simple poles.
"""
struct PartialFractions{T} <: BasisTranslates.Translates{T,T}
    poles   ::  Vector{T}
    support ::  Domain{T}
end

PartialFractions(poles::Vector{T}) where {T} =
    PartialFractions{T}(poles, UnitInterval{T}())

poles(Φ::PartialFractions) = Φ.poles

BasisTranslates.translates_grid(Φ::PartialFractions) = Φ.poles

BasisFunctions.support(Φ::PartialFractions) = Φ.support

kernel_support(Φ::PartialFractions) = support(Φ)
kernel_eval(Φ::PartialFractions, x) = 1/x

function BasisFunctions.dict_norm(Φ::PartialFractions, idx, p::Real)
    @boundscheck checkbounds(Φ, idx)
    @assert (p == 1) || (p==2)
    @assert support(Φ) == 0..1
    c = -Φ.poles[idx]
    if p == 1
        abs(1/c)
    else
        sqrt(-1/(c+1) + 1/c)
    end
end
