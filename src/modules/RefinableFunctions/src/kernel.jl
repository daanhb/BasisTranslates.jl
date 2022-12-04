
"""
Supertype of refinable functions.

A refinable function `φ` satisfies a two-scale relation of the form
`φ(x) = sqrt(2) \\sum_k h[k] φ(2x-k)`
for some set of coefficients `h[k]`.
"""
abstract type Refinable{T} <: BasisTranslates.Kernel end
# type parameter T is supposed to be the eltype of the coefficients

"Supertype of refinable functions with compact support."
abstract type CompactRefinable{T} <: Refinable{T} end

function kernel_support(φ::CompactRefinable{T}) where {T}
    I = support(coefficients(φ))
    T(first(I))..T(last(I))
end

"A refinable function with compact support given by its coefficients."
struct GenericRefinable{T} <: CompactRefinable{T}
    coefficients    ::  VectorSequence{T}
    fun

    GenericRefinable{T}(coefficients) where {T} = new(coefficients)
end

GenericRefinable(coef::AbstractVector) =
    GenericRefinable(VectorSequence(coef))

coefficients(φ::GenericRefinable) = φ.coefficients


## Convenience constructors

Refinable(coef::AbstractVector) = GenericRefinable(coef)
Refinable(coef::VectorSequence) = GenericRefinable(coef)
