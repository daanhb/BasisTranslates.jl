
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

"A refinable function with compact support given by its coefficients."
struct GenericRefinable{T} <: CompactRefinable{T}
    coefficients    ::  VectorSequence{T}
end

GenericRefinable(coef::AbstractVector) =
    GenericRefinable(VectorSequence(coef))

coefficients(φ::GenericRefinable) = φ.coefficients


"The scaling function of the Haar wavelet is the block function on `[0,1]`."
struct HaarScalingFunction{T} <: CompactRefinable{T}
end

HaarScalingFunction() = HaarScalingFunction{Float64}()

coefficients(φ::HaarScalingFunction{T}) where {T} =
    VectorSequence(1/sqrt(T(2))*[1,1])


## Convenience constructors

Refinable(coef::AbstractVector) = GenericRefinable(coef)
Refinable(coef::VectorSequence) = GenericRefinable(coef)
