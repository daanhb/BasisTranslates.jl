
"Supertype of sequences, i.e., bi-infinite arrays of numbers."
abstract type Sequence{T} end

Base.eltype(::Type{<:Sequence{T}}) where {T} = T

"Supertype of compact sequences, with a finite number of non-zero entries."
abstract type CompactSequence{T} <: Sequence{T} end

iscompact(s::Sequence) = false
iscompact(s::CompactSequence) = true

Base.sum(s::CompactSequence) = sum(s[i] for i in support(s))

datavector(s::CompactSequence) = [s[i] for i in support(s)]

"""
A compact vector sequence has a finite number of non-zero entries stored
as a vector.
"""
struct VectorSequence{T} <: CompactSequence{T}
    coefficients    ::  Vector{T}
    I               ::  UnitRange{Int}

    function VectorSequence{T}(coefficients, I) where {T}
        @assert length(I) == length(coefficients)
        new(coefficients, I)
    end
end

# Make a causal sequence by default
VectorSequence(coef::AbstractVector{T}, I = 0:length(coef)-1) where {T} =
    VectorSequence{T}(coef, I)

Base.getindex(s::VectorSequence, i::Int) =
    _getindex(s, i, s.coefficients, s.I)
_getindex(s::VectorSequence{T}, i::Int, coef, I) where {T} =
    i ∈ I ? s.coefficients[i-first(I)+1] : zero(T)

support(s::VectorSequence) = s.I

datavector(s::VectorSequence) = copy(s.coefficients)
unsafe_datavector(s::VectorSequence) = s.coefficients

Base.sum(s::VectorSequence) = sum(unsafe_datavector(s))


"Flip the unit range with respect to 0."
function flipunitrange(I::UnitRange)
    i1, i2 = first(I), last(I)
    (-i2):(-i1)
end

# the reverse of a compact sequence is also a compact sequence
function Base.reverse(s::CompactSequence)
    I = support(s)
    I2 = flipunitrange(I)
    VectorSequence([s[i] for i in reverse(I)], I2)
end

# the adjoint of a sequence is the conjugate of the reversed sequence
function Base.adjoint(s::CompactSequence)
    I = support(s)
    I2 = flipunitrange(I)
    VectorSequence([adjoint(s[i]) for i in reverse(I)], I2)
end
