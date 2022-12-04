
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

function Base.adjoint(A::CompactCirculant)
    seq = conj(reverse(A.seq))
    n = size(A,1)
    CompactCirculant(seq, n, offset = -length(seq)-A.offset+3)
end

# TODO: efficient multiplication

"""
A fat rectangular compact circulant matrix is a subset of a compact
circulant matrix that selects a subset of rows.
"""
struct FatRectangularCirculant{T} <: AbstractArray{T,2}
    C       ::  CompactCirculant{T}
    step    ::  Int
    offset  ::  Int

    function FatRectangularCirculant{T}(C::CompactCirculant; step::Int, offset::Int = 1) where {T}
        @assert 1 <= offset <= step
        new(C, step, offset)
    end
end

FatRectangularCirculant(C::CompactCirculant{T}; step::Int, offset::Int=1) where {T} =
    FatRectangularCirculant{T}(C; step, offset)

function Base.size(A::FatRectangularCirculant)
    m,n = size(A.C)
    div(m, A.step), n
end

Base.getindex(A::FatRectangularCirculant, i::Int, j::Int) =
    A.C[(i-1)*A.step+A.offset,j]

function LinearAlgebra.factorize(A::FatRectangularCirculant)
    C = A.C
    R = StridedRows(size(A,2), step=A.step, offset=A.offset)
    R, C
end

Base.:*(R::BasisTranslates.StridedRows, C::CompactCirculant) =
    FatRectangularCirculant(C; step=step(R), offset=BasisTranslates.offset(R))


"""
A tall rectangular compact circulant matrix is a subset of a compact
circulant matrix that selects a subset of columns.
"""
struct TallRectangularCirculant{T} <: AbstractArray{T,2}
    C       ::  CompactCirculant{T}
    step    ::  Int
    offset  ::  Int

    function TallRectangularCirculant{T}(C::CompactCirculant; step::Int, offset::Int = 1) where {T}
        @assert 1 <= offset <= step
        new(C, step, offset)
    end
end

TallRectangularCirculant(C::CompactCirculant{T}; step::Int, offset::Int=1) where {T} =
    TallRectangularCirculant{T}(C; step, offset)

function Base.size(A::TallRectangularCirculant)
    m,n = size(A.C)
    m, div(n, A.step)
end

Base.getindex(A::TallRectangularCirculant, i::Int, j::Int) =
    A.C[i,(j-1)*A.step+A.offset]

Base.:*(C::CompactCirculant, R::StridedColumns) =
    TallRectangularCirculant(C; step=step(R), offset=BasisTranslates.offset(R))

function LinearAlgebra.factorize(A::TallRectangularCirculant)
    C = A.C
    R = StridedRows(size(A,1), step=A.step, offset=A.offset)
    C, R'
end

function Base.adjoint(A::FatRectangularCirculant)
    R, C = factorize(A)
    C'*R'
end
function Base.adjoint(A::TallRectangularCirculant)
    C, E = factorize(A)
    E'*C'
end

# TODO: efficient multiplication
