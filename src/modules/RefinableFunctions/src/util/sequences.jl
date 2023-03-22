
"Supertype of sequences, i.e., bi-infinite arrays of numbers."
abstract type Sequence{T} end

Base.eltype(::Type{<:Sequence{T}}) where {T} = T
Base.getindex(s::Sequence, I::AbstractRange) = [s[i] for i in I]
Base.getindex(s::Sequence, I::Colon) = s

"Supertype of compact sequences, with a finite number of non-zero entries."
abstract type CompactSequence{T} <: Sequence{T} end

iscompact(s::Sequence) = false
iscompact(s::CompactSequence) = true

Base.sum(s::CompactSequence) = sum(s[i] for i in support(s))

datavector(s::CompactSequence) = [s[i] for i in support(s)]
datalength(s::CompactSequence) = length(support(s))

## Arithmetic with sequences

Base.:*(a::Number, s::CompactSequence) =
    CompactSequence(a*datavector(s), support(s))
Base.:*(s::CompactSequence, a::Number) = a*s
Base.:/(s::CompactSequence, a::Number) = (1/a)*s

function Base.:+(s1::CompactSequence, s2::CompactSequence)
    I1 = support(s1)
    I2 = support(s2)
    I = min(first(I1),first(I2)):max(last(I1),last(I2))
    CompactSequence([s1[i]+s2[i] for i in I], I)
end


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

# Convenience constructor for the abstract type
CompactSequence(coef::AbstractVector, args...) = VectorSequence(coef, args...)

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

"Compute the convolution of two (short) vectors by explicit summation."
function convolve(a::AbstractVector, b::AbstractVector)
    La = length(a)
    Lb = length(b)
    if La < Lb
        return convolve(b, a)
    end
    L = La+Lb-1
    d = zeros(promote_type(eltype(a),eltype(b)), L)
    for i in 1:La
        k1 = max(1,i+1-Lb)
        k2 = min(La,i)
        d[i] = sum(a[k]*b[i+1-k] for k in k1:k2)
    end
    for i in La+1:L
        k1 = i+1-Lb
        k2 = La
        if k1 <= k2     # avoid reduction over empty collection
            d[i] = sum(a[k]*b[i+1-k] for k in k1:k2)
        end
    end
    d
end

"Compute the convolution of two compact sequences."
function convolve(s1::CompactSequence, s2::CompactSequence)
    d = convolve(datavector(s1), datavector(s2))
    i1 = first(support(s1))
    j1 = first(support(s2))
    VectorSequence(d, i1+j1:i1+j1+length(d)-1)
end

"Shift a sequence in time."
shift(s::CompactSequence, i::Int) = CompactSequence(datavector(s), support(s) .+ i)


"The Dirac sequence is zero everywhere except at time `k`, where it is `1`."
struct DiracSequence{T} <: CompactSequence{T}
    k   ::  Int
end

DiracSequence(k::Int = 0) = DiracSequence{Float64}(k)

iscompact(s::DiracSequence) = true
Base.getindex(s::DiracSequence{T}, i::Int) where T =
    i == s.k ? one(T) : zero(T)
support(s::DiracSequence) = s.k:s.k

shift(s::DiracSequence{T}, i::Int) where T = DiracSequence{T}(s.k+i)
