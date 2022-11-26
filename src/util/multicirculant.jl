
"Definition of a block circulant matrix."
const BlockCirculant{T} =
    BlockMatrix{T, Matrix{Circulant{T}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

"Definition of a block diagonal matrix with complex entries."
const ComplexBlockDiagonal{T} =
    BlockMatrix{Complex{T}, Matrix{Diagonal{Complex{T}, Vector{Complex{T}}}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

"Restriction operator of a multi-row circulant matrix."
const RestrictionBlockArray =
    BlockMatrix{Bool, Matrix{U}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}} where {U<:RestrictionArray}

const ColumnBlockArray = Union{
    BlockCirculant{T} where T,
    ComplexBlockDiagonal{T} where T,
    RestrictionBlockArray
}

function Base.:*(A::ColumnBlockArray, x::AbstractVector)
    y = similar(x, promote_type(eltype(A),eltype(x)), axes(A)[1])
    mul!(y, A, x)
end
function Base.:*(A::ColumnBlockArray, x::AbstractBlockVector)
    y = similar(x, promote_type(eltype(A),eltype(x)), axes(A)[1])
    mul!(y, A, x)
end
function Base.:*(A::Adjoint{T,U}, x::AbstractVector) where {T,U<:ColumnBlockArray}
    y = similar(x, promote_type(eltype(A),eltype(x)), Base.OneTo(size(A,1)))
    mul!(y, A, x)
end
function Base.:*(A::Adjoint{T,U}, x::AbstractBlockVector) where {T,U<:ColumnBlockArray}
    y = similar(x, promote_type(eltype(A),eltype(x)), Base.OneTo(size(A,1)))
    mul!(y, A, x)
end

function LinearAlgebra.mul!(y::AbstractBlockVector, A::ColumnBlockArray, x::AbstractVector, α::Number, β::Number)
    @assert α == 1
    @assert iszero(β)
    mul!(y, A, x)
end

function LinearAlgebra.mul!(y::AbstractVector, A::Adjoint{T,U}, x::AbstractBlockVector, α::Number, β::Number) where {T,U<:ColumnBlockArray}
    @assert α == 1
    @assert iszero(β)
    mul!(y, A, x)
end

function LinearAlgebra.mul!(y::AbstractBlockVector, A::ColumnBlockArray, x::AbstractVector)
    m, n = size(A)
    if axes(y)[1] != axes(A)[1]
        throw(DimensionMismatch("block structure of A does not match that of y"))
    end
    if length(x) != n
        throw(DimensionMismatch("second dimension of A, $(n), does not match length of x, $(length(x))"))
    end
    s = blocksize(A, 1)
    for k in 1:s
        mul!(viewblock(y, Block(k)), A[Block(k)], x)
    end
    y
end

function LinearAlgebra.mul!(y::AbstractVector, A::Adjoint{T,U}, x::AbstractBlockVector) where {T,U<:ColumnBlockArray}
    m, n = size(A)
    if length(y) != m
        throw(DimensionMismatch("first dimension of A, $(m), does not match length of y, $(length(y))"))
    end
    if axes(x)[1] != axes(A)[2]
        throw(DimensionMismatch("block structure of A does not match that of x"))
    end
    AA = parent(A)
    s = blocksize(AA, 1)
    tmp = similar(y)
    y .= 0
    for k in 1:s
        mul!(tmp, adjoint(AA[Block(k,1)]), viewblock(x, Block(k)))
        y .+= tmp
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
"""
struct MultiRowCirculantFactorization{T} <: LinearAlgebra.Factorization{T}
    A   ::  MultiRowCirculant{T}
    B   ::  BlockCirculant{T}
    D   ::  ComplexBlockDiagonal{T}
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
