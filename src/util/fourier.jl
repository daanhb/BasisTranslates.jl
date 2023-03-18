############################
# The normalized DFT matrix
############################

using BlockArrays: viewblock

normalized_dft_getindex(n, ::Type{T}, i, j) where {T} =
    exp(2*T(π)*im*(i-1)*(j-1)/n) / sqrt(T(n))

function normalized_dft!(n, y, x, plan!::Nothing = nothing)
    T = eltype(x)
    y[:] = ifft(x)*sqrt(T(n))
    y
end
function normalized_dft!(n, y, x, plan!)
    y[:] = x
    plan! * y
end

function normalized_idft!(n, y, x, iplan!::Nothing = nothing)
    T = eltype(x)
    y[:] = fft(x)/sqrt(T(n))
    y
end
function normalized_idft!(n, y, x, iplan!)
    y[:] = x
    iplan! * y
end


normalized_dft(x::AbstractVector) =
    normalized_dft!(length(x), complex(similar(x)), x)
normalized_idft(x::AbstractVector) =
    normalized_idft!(length(x), complex(similar(x)), x)


"""
The Matrix representing a normalized discrete Fourier transform of
size `n x n`, defined by 'F_{kl} = exp(2πi(k-1)(l-1)/n)/sqrt(n)'.
"""
struct NormalizedDFT{T,P1,P2} <: AbstractArray{Complex{T},2}
    n       ::  Int
    plan!   ::  P1
    iplan!  ::  P2
end
NormalizedDFT(n::Int) = NormalizedDFT{Float64}(n)

function NormalizedDFT{T}(n) where {T <: FFTW.fftwNumber}
    plan! = AbstractFFTs.ScaledPlan(plan_ifft!(zeros(complex(T),n)), sqrt(T(n)))
    iplan! = AbstractFFTs.ScaledPlan(plan_fft!(zeros(complex(T),n)), 1/sqrt(T(n)))
    NormalizedDFT{T}(n, plan!, iplan!)
end

NormalizedDFT{T}(n) where T = NormalizedDFT{T}(n, nothing, nothing)

NormalizedDFT{T}(n, plan!, iplan!) where T =
    NormalizedDFT{T,typeof(plan!),typeof(iplan!)}(n, plan!, iplan!)

const NormalizedDFTAdj{T,P1,P2} = Adjoint{Complex{T},NormalizedDFT{T,P1,P2}}

Base.size(A::NormalizedDFT) = (A.n, A.n)

@inline function Base.getindex(A::NormalizedDFT{T}, k::Int, l::Int) where {T}
    @boundscheck checkbounds(A, k, l)
    normalized_dft_getindex(A.n, T, k, l)
end

function check_dimensions(m, n, y, x)
    if length(y) != m
        throw(DimensionMismatch("first dimension of A, $(m), does not match length of y, $(length(y))"))
    end
    if length(x) != n
        throw(DimensionMismatch("second dimension of A, $(n), does not match length of x, $(length(x))"))
    end
end

function LinearAlgebra.mul!(y::AbstractVector, F::NormalizedDFT, x::AbstractVector, α::Number, β::Number)
    mul!(y, F, x)
    if α != 1
        y .= α .* y
    end
    if !iszero(β)
        y .= y .+ β .* x
    end
    y
end
function LinearAlgebra.mul!(y::AbstractVector, F::NormalizedDFT, x::AbstractVector)
    m, n = size(F)
    check_dimensions(m, n, y, x)
    normalized_dft!(n, y, x, F.plan!)
end

function LinearAlgebra.mul!(y::AbstractVector, A::NormalizedDFTAdj, x::AbstractVector, α::Number, β::Number)
    mul!(y, A, x)
    if α != 1
        y .= α .* y
    end
    if !iszero(β)
        y .= y .+ β .* x
    end
    y
end
function LinearAlgebra.mul!(y::AbstractVector, A::NormalizedDFTAdj, x::AbstractVector)
    n = size(A, 1)
    normalized_idft!(n, y, x, parent(A).iplan!)
end

LinearAlgebra.inv(A::NormalizedDFT) = adjoint(A)
LinearAlgebra.pinv(A::NormalizedDFT) = inv(A)

LinearAlgebra.inv(A::NormalizedDFTAdj) = adjoint(A)
LinearAlgebra.pinv(A::NormalizedDFTAdj) = inv(A)


##################################
# A block-diagonal Fourier matrix
##################################

"""
A block matrix with `s x s` blocks, each of which has size `n x n`.
The diagonal blocks of this matrix are normalized DFT matrices.
"""
struct BlockFourierMatrix{T,A} <: AbstractBlockArray{Complex{T},2}
    s       ::  Int
    n       ::  Int
    block   ::  A
    axes1   ::  BlockedUnitRange{Vector{Int64}}

    function BlockFourierMatrix{T,A}(s::Int, n::Int, block::NormalizedDFT) where {T,A}
        axes1 = blockedrange([n for i in 1:s])
        new(s, n, block, axes1)
    end
end

BlockFourierMatrix{T}(s::Int, n::Int) where T =
    BlockFourierMatrix{T}(s, n, NormalizedDFT{T}(n))
BlockFourierMatrix{T}(s::Int, n::Int, block::A) where {T,A} =
    BlockFourierMatrix{T,A}(s, n, block)

const BlockFourierMatrixAdj{T,A} = Adjoint{Complex{T},BlockFourierMatrix{T,A}}

Base.axes(A::BlockFourierMatrix) = (A.axes1, A.axes1)

function Base.getindex(A::BlockFourierMatrix, i::Int, j::Int)
    @boundscheck checkbounds(A, i, j)
    I1,I2 = findblockindex.(axes(A), (i,j))
    if I1.I == I2.I
        k = I1.α[1]
        l = I2.α[1]
        A.block[k,l]
    else
        zero(eltype(A))
    end
end

function Base.:*(F::BlockFourierMatrix, x::AbstractBlockVector)
    y = similar(x, promote_type(eltype(F),eltype(x)))
    LinearAlgebra.mul!(y, F, x)
end
function Base.:*(F::BlockFourierMatrix, x::AbstractVector)
    x2 = PseudoBlockVector(x, blocksizes(F)[2])
    y = similar(x2, promote_type(eltype(F),eltype(x2)))
    LinearAlgebra.mul!(y, F, x2)
end

# five argument call: catch and simply (throw an error if not applicable)
function LinearAlgebra.mul!(y::AbstractVector, F::BlockFourierMatrix, x::AbstractVector, α::Number, β::Number)
    @assert (α == 1) && iszero(β)
    mul!(y, F, x)
end
# convert regular arrays to block arrays
function LinearAlgebra.mul!(y::AbstractVector, F::BlockFourierMatrix, x::AbstractVector)
    x2 = PseudoBlockVector(x, blocksizes(F)[2])
    y2 = PseudoBlockVector(y, blocksizes(F)[1])
    mul!(y2, F, x2)
    y
end
function LinearAlgebra.mul!(y::AbstractVector, F::BlockFourierMatrix, x::AbstractBlockVector)
    y2 = PseudoBlockVector(y, blocksizes(F)[1])
    mul!(y2, F, x)
    y
end
function LinearAlgebra.mul!(y::AbstractBlockVector, F::BlockFourierMatrix, x::AbstractVector)
    x2 = PseudoBlockVector(x, blocksizes(F)[2])
    mul!(y, F, x2)
end

# and now comes the real implementation
function LinearAlgebra.mul!(y::AbstractBlockVector, F::BlockFourierMatrix, x::AbstractBlockVector)
    @assert axes(x)[1] == axes(y)[1] == axes(F)[2]
    s, n = F.s, F.n
    for k in 1:s
        mul!(viewblock(y, Block(k)), F.block, viewblock(x, Block(k)))
    end
    y
end

function Base.:*(A::BlockFourierMatrixAdj, x::AbstractBlockVector)
    y = similar(x, promote_type(eltype(A),eltype(x)))
    LinearAlgebra.mul!(y, A, x)
end
function Base.:*(A::BlockFourierMatrixAdj, x::AbstractVector)
    F = parent(A)
    x2 = PseudoBlockVector(x, blocksizes(F)[1])
    y = similar(x2, promote_type(eltype(A),eltype(x2)))
    LinearAlgebra.mul!(y, A, x2)
end

# catch the five-argument call
function LinearAlgebra.mul!(y::AbstractVector, A::BlockFourierMatrixAdj, x::AbstractVector, α::Number, β::Number)
    @assert (α == 1) && iszero(β)
    mul!(y, A, x)
end

# convert vectors to block vectors if necessary
function LinearAlgebra.mul!(y::AbstractVector, A::BlockFourierMatrixAdj, x::AbstractVector)
    F = parent(A)
    x2 = PseudoBlockVector(x, blocksizes(F)[1])
    y2 = PseudoBlockVector(y, blocksizes(F)[2])
    mul!(y2, A, x2)
    y
end
function LinearAlgebra.mul!(y::AbstractBlockVector, A::BlockFourierMatrixAdj, x::AbstractVector)
    F = parent(A)
    x2 = PseudoBlockVector(x, blocksizes(F)[1])
    mul!(y, A, x2)
end
function LinearAlgebra.mul!(y::AbstractVector, A::BlockFourierMatrixAdj, x::AbstractBlockVector)
    F = parent(A)
    y2 = PseudoBlockVector(y, blocksizes(F)[2])
    mul!(y2, A, x)
    y
end

function LinearAlgebra.mul!(y::AbstractBlockVector, A::BlockFourierMatrixAdj, x::AbstractBlockVector)
    F = parent(A)
    @assert axes(x)[1] == axes(y)[1] == axes(F)[2]
    s, n = F.s, F.n
    for k in 1:s
        mul!(viewblock(y, Block(k)), parent(A).block', viewblock(x, Block(k)))
    end
    y
end

LinearAlgebra.inv(A::BlockFourierMatrix) = adjoint(A)
LinearAlgebra.pinv(A::BlockFourierMatrix) = inv(A)

LinearAlgebra.inv(A::BlockFourierMatrixAdj) = adjoint(A)
LinearAlgebra.pinv(A::BlockFourierMatrixAdj) = inv(A)
