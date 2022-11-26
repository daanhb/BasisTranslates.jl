"Definition of a column-block circulant matrix."
const BlockCirculant{T} =
    BlockMatrix{T, Matrix{Circulant{T}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

"Definition of a column-block diagonal matrix."
const BlockDiagonal{T} =
    BlockMatrix{T, Matrix{Diagonal{T, Vector{T}}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

"Block row selection operator of a multi-row circulant matrix."
const RestrictionBlockArray =
    BlockMatrix{Bool, Matrix{U}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}} where {U<:RestrictionArray}

const ColumnBlockArray = Union{
    BlockCirculant{T} where T,
    BlockDiagonal{T} where T,
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
