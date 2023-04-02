
"""
A multi-circulant matrix of intertwined circulant matrices, row by row.
"""
struct MultiRowCirculant{T} <: AbstractMatrix{T}
    arrays  ::  Vector{Circulant{T}}
end

nblocks(A::MultiRowCirculant) = length(A.arrays)
blockdim(A::MultiRowCirculant) = size(A.arrays[1])[1]

function Base.size(A::MultiRowCirculant)
    s, n = nblocks(A), blockdim(A)
    n*s,n
end

function Base.getindex(A::MultiRowCirculant, i, j)
    s = nblocks(A)
    k, l = divrem(i-1, s)
    A.arrays[l+1][k+1,j]
end

Base.copy(A::MultiRowCirculant) = MultiRowCirculant([copy(a) for a in A.arrays])

"Convert the given matrix into a block-column circulant matrix."
function blockcirculant(A::MultiRowCirculant)
    s = nblocks(A)
    mortar(reshape(A.arrays, s, 1))
end

"Permutation to convert the given matrix into a block-column circulant matrix."
function rowpermutation(A::MultiRowCirculant{T}) where T
    s, n = nblocks(A), blockdim(A)
    ops = [StridedRows{Complex{T}}(s*n; step = s, offset = i) for i in 1:s]
    mortar(reshape(ops, s, 1))
end

function LinearAlgebra.mul!(y::AbstractVector, A::MultiRowCirculant, x::AbstractVector)
    s, n = nblocks(A), blockdim(A)
    for i in 1:s
        yk = A.arrays[i]*x
        y[i:s:end] = yk
    end
    y
end

const RowPermutationArray{T} = BlockMatrix{T, Matrix{StridedRows{T}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

"""
Factorization of a multi-row circulant matrix.

The multi-row ciculant matrix `A` corresponds to a column-block circulant
matrix `B = Π*A` and `A = Π'*B`.
The matrix `B` is block-diagonalized using block Fourier, `P'*B*F = D`.

The complete factorization satisfies `A = Π'*P*D*F'`.
"""
struct MultiRowCirculantFactorization{T} <: LinearAlgebra.Factorization{Complex{T}}
    A       ::  MultiRowCirculant{T}
    B       ::  BlockCirculant{T}
    D       ::  BlockDiagonal{Complex{T}}
    Dpinv   ::  BlockDiagonalAdj{Complex{T}}
    F       ::  NormalizedDFT{T}
    P       ::  BlockFourierMatrix{T}
    Π       ::  RowPermutationArray{Complex{T}}
end

nblocks(F::MultiRowCirculantFactorization) = nblocks(F.A)

function LinearAlgebra.factorize(A::MultiRowCirculant{T}) where {T}
    s, n = nblocks(A), blockdim(A)
    B = blockcirculant(A)
    diagonals = [Diagonal(eigvals(bl)) for bl in eachblock(B)]
    D = mortar(reshape(diagonals, s, 1))
    Dpinv = pinv(D)
    F = NormalizedDFT{T}(n)
    P = BlockFourierMatrix{T}(s, n)
    Π = rowpermutation(A)
    MultiRowCirculantFactorization(A, B, D, Dpinv, F, P, Π)
end

import LinearAlgebra: \

# solve x = A \ b
function (\)(Fact::MultiRowCirculantFactorization, b::AbstractVector)
    c = Fact.P' * (Fact.Π * b)
    y = Fact.Dpinv * c
    x = Fact.F * y
    isreal(b) ? real(x) : x
end
