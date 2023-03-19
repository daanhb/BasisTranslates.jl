
"""
Supertype of refinable functions.

A refinable function `φ` satisfies a two-scale relation of the form
`φ(x) = sqrt(2) \\sum_k h[k] φ(2x-k)`
for some set of coefficients `h[k]`.
"""
abstract type Refinable{T} <: Kernel end
# type parameter T is supposed to be the eltype of the coefficients

numtype(::Type{<:Refinable{T}}) where T = T

isrefinable(φ::Refinable) = true
isrefinable(φ::Kernel) = false

"Return the coefficients of the two-scale relation as a sequence."
refinable_coeff(φ::BasisTranslates.PeriodizedKernel) =
    refinable_coeff(parent_kernel(φ))

"""
A refinable function is only defined up to a normalization. The normalization
can be fixed by the value of the integral of the function (its first moment).

Since refinable functions form a partition of unity, the value of the moment
also equals the sum of the values of the function at equispaced points with
grid spacing 1.
"""
refinable_moment(φ::Refinable) = one(numtype(φ))


"Supertype of refinable functions with compact support."
abstract type CompactRefinable{T} <: Refinable{T} end

function kernel_support(φ::CompactRefinable{T}) where {T}
    I = support(refinable_coeff(φ))
    T(first(I))..T(last(I))
end

"""
A refinable function with compact support given by its coefficients.
The function is normalized to have a given integral.
"""
struct GenericRefinable{T} <: CompactRefinable{T}
    coefficients    ::  VectorSequence{T}
    moment          ::  T
end

GenericRefinable(coef::AbstractVector, args...) =
    GenericRefinable(VectorSequence(coef), args...)
GenericRefinable(coefficients::VectorSequence{T}, moment = one(T)) where T =
    GenericRefinable{T}(coefficients, moment)

"Return the coefficients of the refinement equation of the given refinable function."
refinable_coeff(φ::GenericRefinable) = φ.coefficients
"Return the first moment of the refinable function."
refinable_moment(φ::GenericRefinable) = φ.moment

"Alias for `refinable_coeff`."
refinable_coefficients(φ) = refinable_coeff(φ)

"""
A generic refinable function whose evaluation is computed in a fine
grid upon construction and stored for quick evaluation later on.
"""
struct EvaluatedRefinable{T} <: CompactRefinable{T}
    parent  ::  Kernel
    fun
    function EvaluatedRefinable{T}(φ::Kernel; levels = 6) where {T}
        @assert isrefinable(φ)
        t, vals = eval_dyadic(refinable_coeff(φ), levels, refinable_moment(φ))
        fun = Expansion(RegularBSplines(t, degree=1), vals)
        new(φ, fun)
    end
end

EvaluatedRefinable(φ::Refinable{T}; options...) where {T} =
    EvaluatedRefinable{T}(φ; options...)
EvaluatedRefinable(φ::Kernel; options...) =
    EvaluatedRefinable{eltype(refinable_coeff(φ))}(φ; options...)

Base.parent(φ::EvaluatedRefinable) = φ.parent

refinable_coeff(φ::EvaluatedRefinable) = refinable_coeff(parent(φ))

kernel_eval(φ::EvaluatedRefinable, x) = _kernel_eval(φ, x, φ.fun)
_kernel_eval(φ::EvaluatedRefinable, x, fun) = fun(x)

## Convenience constructors

Refinable(coef::AbstractVector) = GenericRefinable(coef)
Refinable(coef::VectorSequence) = GenericRefinable(coef)
