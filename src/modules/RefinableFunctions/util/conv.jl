
"Periodize a vector by summing over segments of length `n`."
function periodize_sequence(seq::AbstractVector, n::Int)
    if length(seq) > n
        periodic_seq = similar(seq, n)
        for i in 1:n
            periodic_seq[i] = sum(seq[i:n:length(seq)])
        end
        periodic_seq
    else
        seq
    end
end

"""
A compact circulant array is a circulant array that is also bandlimited,
in the sense that all its rows have only `L` nonzero entries.
The array is defined by its dimensions, a vector of length `L` and an offset
for the position of that vector on the first row.
"""
struct CompactCirculant{T} <: AbstractArray{T,2}
    seq     ::  Vector{T}
    n       ::  Int
    offset  ::  Int

    function CompactCirculant{T}(seq::AbstractVector, n; offset = 1) where {T}
        @assert length(seq) <= n
        @assert 1 <= offset <= n
        new(seq, n, offset)
    end
end

CompactCirculant(seq::AbstractVector{T}, n::Int; offset = 1) where {T} =
    CompactCirculant{T}(seq, n; offset)

Base.size(A::CompactCirculant) = (A.n, A.n)

function Base.getindex(A::CompactCirculant{T}, i::Int, j::Int) where {T}
    @boundscheck checkbounds(A, i, j)
    idx = mod(j-i-A.offset+1,A.n) + 1
    1 <= idx <= length(A.seq) ? A.seq[idx] : zero(T)
end

# TODO: efficient multiplication

"""
A rectangular compact circulant matrix is a subset of a compact
circulant matrix.
"""
struct RectangularCirculant{T} <: AbstractArray{T,2}
    A       ::  CompactCirculant{T}
    step    ::  Int
    offset  ::  Int

    function RectangularCirculant{T}(A::CompactCirculant, step::Int; offset::Int = 1) where {T}
        @assert 1 <= offset <= step
        new(A, step, offset)
    end
end

RectangularCirculant(A::CompactCirculant{T}, step::Int; offset::Int=1) where {T} =
    RectangularCirculant{T}(A, step; offset)

function Base.size(A::RectangularCirculant)
    m,n = size(A.A)
    div(m, A.step), n
end

Base.getindex(A::RectangularCirculant, i::Int, j::Int) =
    A.A[(i-1)*A.step+A.offset,j]

Base.:*(R::BasisTranslates.StridedRows, A::CompactCirculant) =
    RectangularCirculant(A, step(R); offset=BasisTranslates.offset(R))

# TODO: efficient multiplication
