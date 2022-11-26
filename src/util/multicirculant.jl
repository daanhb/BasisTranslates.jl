
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

"Convert the given matrix into a block circulant matrix."
function blockcirculant(A::MultiRowCirculant)
    s = nblocks(A)
    mortar(reshape(A.arrays, s, 1))
end

"Convert the given matrix into a block circulant matrix."
function blockrowselection(A::MultiRowCirculant)
    s, n = nblocks(A), blockdim(A)
    L = s*n
    shifted_grids = [k:s:L for k in 1:s]
    ops = StridedRows.(L, shifted_grids)
    mortar(reshape(ops, s, 1))
end


const MultiRowPermutationArray = BlockMatrix{Bool, Matrix{StridedRows}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

"""
Factorization of a multi-row circulant matrix.

The multi-row ciculant matrix `A` corresponds to a column-block circulant
matrix `B = R*A` and `A = R'*B`.
The matrix `B` is block-diagonalized using block Fourier, `P'*B*F = D`.

The complete factorization satisfies `A = R'*P*D*F`.
"""
struct MultiRowCirculantFactorization{T} <: LinearAlgebra.Factorization{T}
    A   ::  MultiRowCirculant{T}
    B   ::  BlockCirculant{T}
    D   ::  BlockDiagonal{Complex{T}}
    F   ::  NormalizedDFT{T}
    P   ::  BlockFourierMatrix{T}
    R   ::  MultiRowPermutationArray
end

nblocks(F::MultiRowCirculantFactorization) = nblocks(F.A)

function LinearAlgebra.factorize(A::MultiRowCirculant{T}) where {T}
    s, n = nblocks(A), blockdim(A)
    B = blockcirculant(A)
    diagonals = [Diagonal(eigvals(bl)) for bl in eachblock(B)]
    D = mortar(reshape(diagonals, s, 1))
    F = NormalizedDFT{T}(n)
    P = BlockFourierMatrix{T}(s, n)
    R = blockrowselection(A)
    MultiRowCirculantFactorization(A, B, D, F, P, R)
end
