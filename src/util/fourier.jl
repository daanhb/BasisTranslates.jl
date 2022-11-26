############################
# The normalized DFT matrix
############################

using BlockArrays: viewblock

normalized_dft_getindex(n, ::Type{T}, i, j) where {T} =
    exp(2*T(π)*im*(i-1)*(j-1)/n) / sqrt(T(n))

function normalized_dft!(n, ::Type{T}, y, x) where {T}
    y[:] = ifft(x)*sqrt(T(n))
    y
end

function normalized_idft!(n, ::Type{T}, y, x) where {T}
    y[:] = fft(x)/sqrt(T(n))
    y
end

normalized_dft(x::AbstractVector) =
    normalized_dft!(length(x), eltype(x), complex(similar(x)), x)
normalized_idft(x::AbstractVector) =
    normalized_idft!(length(x), eltype(x), complex(similar(x)), x)


"""
The Matrix representing a normalized discrete Fourier transform of
size `n x n`, defined by 'F_{kl} = exp(2πi(k-1)(l-1)/n)/sqrt(n)'.
"""
struct NormalizedDFT{T} <: AbstractArray{Complex{T},2}
    n       ::  Int
end
NormalizedDFT(n::Int) = NormalizedDFT{Float64}(n)

const NormalizedDFTAdj{T} = Adjoint{Complex{T},NormalizedDFT{T}}

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

function LinearAlgebra.mul!(y::AbstractVector, F::NormalizedDFT{T}, x::AbstractVector, α::Number, β::Number) where {T}
    m, n = size(F)
    check_dimensions(m, n, y, x)
    normalized_dft!(n, T, y, x)
    if α != 1
        y .= α .* y
    end
    if !iszero(β)
        y .= y .+ β .* x
    end
    y
end

function LinearAlgebra.mul!(y::AbstractVector, A::NormalizedDFTAdj{T}, x::AbstractVector, α::Number, β::Number) where {T}
    n = size(A, 1)
    normalized_idft!(n, T, y, x)
    if α != 1
        y .= α .* y
    end
    if !iszero(β)
        y .= y .+ β .* x
    end
    y
end


##################################
# A block-diagonal Fourier matrix
##################################

"""
A block matrix with `s x s` blocks, each of which has size `n x n`.
The diagonal blocks of this matrix are normalized DFT matrices.
"""
struct BlockFourierMatrix{T} <: AbstractBlockArray{Complex{T},2}
    s   ::  Int
    n   ::  Int
    axes1   :: BlockedUnitRange{Vector{Int64}}

    function BlockFourierMatrix{T}(s::Int, n::Int) where {T}
        axes1 = blockedrange([n for i in 1:s])
        new(s, n, axes1)
    end
end

const BlockFourierMatrixAdj{T} = Adjoint{Complex{T},BlockFourierMatrix{T}}

Base.axes(A::BlockFourierMatrix) = (A.axes1, A.axes1)

function Base.getindex(A::BlockFourierMatrix{T}, i::Int, j::Int) where {T}
    @boundscheck checkbounds(A, i, j)
    I1,I2 = findblockindex.(axes(A), (i,j))
    if I1.I == I2.I
        k = I1.α[1]
        l = I2.α[1]
        normalized_dft_getindex(A.n, T, k, l)
    else
        zero(Complex{T})
    end
end

function Base.:*(F::BlockFourierMatrix, x::AbstractBlockVector)
    y = similar(x, promote_type(eltype(F),eltype(x)))
    LinearAlgebra.mul!(y, F, x)
end
function LinearAlgebra.mul!(y::AbstractBlockVector, F::BlockFourierMatrix{T}, x::AbstractBlockVector) where {T}
    @assert axes(x)[1] == axes(y)[1] == axes(F)[2]
    s, n = F.s, F.n
    for k in 1:s
        normalized_dft!(n, T, viewblock(y, Block(k)), viewblock(x, Block(k)))
    end
    y
end

function Base.:*(A::BlockFourierMatrixAdj, x::AbstractBlockVector)
    y = similar(x, promote_type(eltype(A),eltype(x)))
    LinearAlgebra.mul!(y, A, x)
end
function LinearAlgebra.mul!(y::AbstractBlockVector, A::BlockFourierMatrixAdj{T}, x::AbstractBlockVector) where {T}
    F = parent(A)
    @assert axes(x)[1] == axes(y)[1] == axes(F)[2]
    s, n = F.s, F.n
    for k in 1:s
        normalized_idft!(n, T, viewblock(y, Block(k)), viewblock(x, Block(k)))
    end
    y
end
