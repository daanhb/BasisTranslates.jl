
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
end

GenericRefinable(coef::AbstractVector) =
    GenericRefinable(VectorSequence(coef))

coefficients(φ::GenericRefinable) = φ.coefficients


"""
A generic refinable function whose evaluation is computed in a fine
grid upon construction and stored for quick evaluation later on.
"""
struct EvaluatedRefinable{T} <: CompactRefinable{T}
    parent  ::  Refinable{T}
    fun
    function EvaluatedRefinable{T}(φ::Refinable; levels = 6) where {T}
        t, vals = eval_dyadic(coefficients(φ), levels)
        fun = Expansion(RegularBSplines(t, degree=1), vals)
        new(φ, fun)
    end
end

EvaluatedRefinable(φ::Refinable{T}) where {T} = EvaluatedRefinable{T}(φ)

Base.parent(φ::EvaluatedRefinable) = φ.parent

coefficients(φ::EvaluatedRefinable) = coefficients(parent(φ))

kernel_eval(φ::EvaluatedRefinable, x) = _kernel_eval(φ, x, φ.fun)
_kernel_eval(φ::EvaluatedRefinable, x, fun) = fun(x)

## Convenience constructors

Refinable(coef::AbstractVector) = GenericRefinable(coef)
Refinable(coef::VectorSequence) = GenericRefinable(coef)
