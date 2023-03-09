
"Supertype of translates of refinable kernels."
abstract type Refinables{T,EXT} <: BasisTranslates.Translates{T,T,EXT} end

"A basis of translates of a (periodized) refinable function."
struct PeriodicRefinables{T} <: Refinables{T,:unitperiodic}
    refinable   ::  Kernel
    n           ::  Int
end

PeriodicRefinables(φ::Refinable{T}, n::Int) where {T} =
    PeriodicRefinables{T}(φ, n)
PeriodicRefinables(φ::Kernel, n::Int) =
    PeriodicRefinables{eltype(refinable_coeff(φ))}(φ, n)

Base.size(Φ::PeriodicRefinables) = (Φ.n,)

Base.similar(Φ::PeriodicRefinables{T}, ::Type{T}, n::Int) where {T} =
    PeriodicRefinables(Φ.refinable, n)

parent_kernel(Φ::PeriodicRefinables) = Φ.refinable

refinable_coeff(Φ::PeriodicRefinables) = refinable_coeff(parent_kernel(Φ))

function refine(f::Expansion{S,T,B}) where {S,T,B<:PeriodicRefinables}
    n = length(f)
    basis2 = similar(dictionary(f), 2n)
    coef = refinable_coeff(parent_kernel(basis2))
    I = support(coef)
    A = CompactCirculant(datavector(coef), 2n; offset = first(I)+1)
    R = RestrictionArray(2n, 1:2:2n)
    C = R*A
    Expansion(basis2, C'*coefficients(f)*sqrt(S(2)))
end
