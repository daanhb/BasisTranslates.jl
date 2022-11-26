
"Supertype of rectangular arrays that select a subset of a vector."
abstract type RestrictionArray <: AbstractArray{Bool,2} end

const ExtensionArray = Adjoint{Bool,<:RestrictionArray}

Base.size(A::RestrictionArray) = (length(rowselection(A)), columnlength(A))

@inline function Base.getindex(A::RestrictionArray, i::Int, j::Int)
    @boundscheck checkbounds(A, i, j)
    j == rowselection(A)[i]
end

function LinearAlgebra.mul!(y::AbstractVector, R::RestrictionArray, x::AbstractVector, α::Number, β::Number)
    m, n = size(R)
    check_dimensions(m, n, y, x)
    @assert α == 1
    @assert iszero(β)
    for (k,i) in enumerate(rowselection(R))
        y[k] = x[i]
    end
    y
end

function LinearAlgebra.mul!(y::AbstractVector, E::ExtensionArray, x::AbstractVector, α::Number, β::Number)
    m, n = size(E)
    if length(y) != m
        throw(DimensionMismatch("first dimension of A, $(m), does not match length of y, $(length(y))"))
    end
    if length(x) != n
        throw(DimensionMismatch("second dimension of A, $(n), does not match length of x, $(length(x))"))
    end
    @assert α == 1
    @assert iszero(β)
    y .= 0
    for (k,i) in enumerate(rowselection(parent(E)))
        y[i] = x[k]
    end
    y
end


"Selecting a step-range of elements."
struct StridedRows <: RestrictionArray
    len         ::  Int
    selection   ::  StepRange{Int,Int}
end

columnlength(A::StridedRows) = A.len
rowselection(A::StridedRows) = A.selection
