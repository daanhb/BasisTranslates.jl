
"Supertype of rectangular arrays that select a subset of a vector."
abstract type RestrictionArray <: AbstractArray{Bool,2} end

# The adjoint of a restriction is an extension by zero padding
const ExtensionArray = Adjoint{Bool,<:RestrictionArray}

Base.size(A::RestrictionArray) = (length(rowselection(A)), columnlength(A))

@inline function Base.getindex(A::RestrictionArray, i::Int, j::Int)
    @boundscheck checkbounds(A, i, j)
    j == rowselection(A)[i]
end

function LinearAlgebra.mul!(y::AbstractVector, R::RestrictionArray, x::AbstractVector, α::Number, β::Number)
    @assert (α == 1) && iszero(β)
    mul!(y, R, x)
end

function LinearAlgebra.mul!(y::AbstractVector, R::RestrictionArray, x::AbstractVector)
    m, n = size(R)
    check_dimensions(m, n, y, x)
    for (k,i) in enumerate(rowselection(R))
        y[k] = x[i]
    end
    y
end

function LinearAlgebra.mul!(y::AbstractVector, E::ExtensionArray, x::AbstractVector, α::Number, β::Number)
    @assert (α == 1) && iszero(β)
    mul!(y, E, x)
end

function LinearAlgebra.mul!(y::AbstractVector, E::ExtensionArray, x::AbstractVector)
    m, n = size(E)
    check_dimensions(m, n, y, x)
    y .= 0
    for (k,i) in enumerate(rowselection(parent(E)))
        y[i] = x[k]
    end
    y
end

LinearAlgebra.pinv(A::RestrictionArray) = adjoint(A)
LinearAlgebra.pinv(A::ExtensionArray) = adjoint(A)


"Selecting a step-range of elements."
struct StridedRows <: RestrictionArray
    len         ::  Int
    selection   ::  StepRange{Int,Int}
end

# Make a StridedRows selection with range offset:step:len
StridedRows(len::Int, step::Int, offset = 1) = StridedRows(len, offset:step:len)

columnlength(A::StridedRows) = A.len
rowselection(A::StridedRows) = A.selection


"Selecting a unit range of elements."
struct ContiguousRows <: RestrictionArray
    len         ::  Int
    selection   ::  UnitRange{Int}
end

columnlength(A::ContiguousRows) = A.len
rowselection(A::ContiguousRows) = A.selection
