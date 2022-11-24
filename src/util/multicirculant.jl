
"""
A multi-circulant matrix of intertwined circulant matrices, row by row.
"""
struct MultiRowCirculant{T} <: AbstractMatrix{T}
    arrays  ::  Vector{Circulant{T}}
end

nblocks(A::MultiRowCirculant) = length(A.arrays)

function Base.size(A::MultiRowCirculant)
    n = size(A.arrays[1])[1]
    s = nblocks(A)
    n*s,n
end

function Base.getindex(A::MultiRowCirculant, i, j)
    s = nblocks(A)
    k, l = divrem(i-1, s)
    A.arrays[l+1][k+1,j]
end

Base.copy(A::MultiRowCirculant) = MultiRowCirculant([copy(a) for a in A.arrays])

function blockcirculant(A::MultiRowCirculant)
    s = nblocks(A)
    mortar(reshape(A.arrays, s, 1))
end
