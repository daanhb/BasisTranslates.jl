
"Supertype of sequences, i.e., bi-infinite arrays of numbers."
abstract type Sequence{T} end

Base.eltype(::Type{<:Sequence{T}}) where {T} = T

"Supertype of compact sequences, with a finite number of non-zero entries."
abstract type CompactSequence{T} <: Sequence{T} end

iscompact(s::Sequence) = false
iscompact(s::CompactSequence) = true

Base.sum(s::CompactSequence) = sum(s[i] for i ∈ support(s))

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

vector(s::VectorSequence) = copy(s.coefficients)
unsafe_vector(s::VectorSequence) = s.coefficients

Base.sum(s::VectorSequence) = sum(unsafe_vector(s))
