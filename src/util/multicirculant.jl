
"Definition of a block circulant matrix."
const BlockCirculant{T} =
    BlockMatrix{T, Matrix{Circulant{T}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

"Definition of a block diagonal matrix with complex entries."
const ComplexBlockDiagonal{T} =
    BlockMatrix{Complex{T}, Matrix{Diagonal{Complex{T}, Vector{Complex{T}}}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

const ColumnBlockArray{T} = Union{BlockCirculant{T}, ComplexBlockDiagonal{T}}

function LinearAlgebra.mul!(y::AbstractBlockVector, A::ColumnBlockArray, x::AbstractBlockVector, α::Number, β::Number)
    @assert α == 1
    @assert iszero(β)
    mul!(y, A, x)
end

function LinearAlgebra.mul!(y::AbstractBlockVector, A::ColumnBlockArray{T}, x::AbstractBlockVector) where {T}
    s, _ = blocksize(A)
    @assert axes(x)[1] == axes(y)[1] == axes(F)[2]
    for k in 1:s
        mul!(viewblock(y, Block(k)), A[Block(k)], viewblock(x, Block(k)))
    end
    y
end

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



"""
Factorization of a multi-row circulant matrix.
"""
struct MultiRowCirculantFactorization{T} <: LinearAlgebra.Factorization{T}
    A   ::  MultiRowCirculant{T}
    B   ::  BlockCirculant{T}
    D   ::  ComplexBlockDiagonal{T}
    F   ::  NormalizedDFT{T}
    P   ::  BlockFourierMatrix{T}
end

nblocks(F::MultiRowCirculantFactorization) = nblocks(F.A)

function LinearAlgebra.factorize(A::MultiRowCirculant{T}) where {T}
    s, n = nblocks(A), blockdim(A)
    B = blockcirculant(A)
    diagonals = [Diagonal(eigvals(bl)) for bl in eachblock(B)]
    D = mortar(reshape(diagonals, s, 1))
    F = NormalizedDFT{T}(n)
    P = BlockFourierMatrix{T}(s, n)
    MultiRowCirculantFactorization(A, B, D, F, P)
end
