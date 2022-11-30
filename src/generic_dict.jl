
"Dictionary of translates of a given kernel function."
struct KernelTranslates{S,T,K} <: Translates{S,T}
    kernel      ::  K
    centers     ::  Vector{S}
    support     ::  Domain{S}
end

KernelTranslates(kernel, centers::AbstractVector) =
    KernelTranslates(kernel, centers, UnitInterval{eltype(centers)}())
KernelTranslates(kernel, centers::AbstractVector{S}, support::Domain{T}) where {S,T} =
    KernelTranslates{promote_type(S,T)}(kernel, centers, support)
KernelTranslates{S}(kernel, centers, support) where {S} =
    KernelTranslates{S,S}(kernel, centers, support)
KernelTranslates{S,T}(kernel, centers, support) where {S,T} =
    KernelTranslates{S,T,typeof(kernel)}(kernel, centers, support)

kernel(Φ::KernelTranslates) = Φ.kernel

translates_grid(Φ::KernelTranslates) = Φ.centers

support(Φ::KernelTranslates) = Φ.support


"Dictionary of periodized translates of a given kernel function."
struct PeriodicKernelTranslates{S,T,K} <: PeriodicTranslates{S,T}
    kernel      ::  K
    n           ::  Int
    scaling     ::  Bool
end

PeriodicKernelTranslates(kernel, args...) =
    PeriodicKernelTranslates{Float64}(kernel, args...)
PeriodicKernelTranslates{T}(kernel, args...) where {T} =
    PeriodicKernelTranslates{T,T}(kernel, args...)
PeriodicKernelTranslates{S,T}(kernel, args...) where {S,T} =
    PeriodicKernelTranslates{S,T,typeof(kernel)}(kernel, args...)
# By default we assume linear scaling
PeriodicKernelTranslates{S,T,K}(kernel::K, n::Int) where {S,T,K} =
    PeriodicKernelTranslates(kernel, n, true)

kernel(Φ::PeriodicKernelTranslates) = Φ.kernel

Base.size(Φ::PeriodicKernelTranslates) = (Φ.n,)

BasisTranslates.linearscaling(Φ::PeriodicKernelTranslates) = Φ.scaling

kernel_eval(Φ::PeriodicKernelTranslates, x) =
    kernel_eval(kernel(Φ), x)
kernel_eval(kernel, x) = kernel(x)

kernel_eval_derivative(Φ::PeriodicKernelTranslates, order, x) =
    kernel_eval_derivative(kernel(Φ), order, x)

kernel_support(Φ::PeriodicKernelTranslates) = kernel_support(kernel(Φ))

kernel_support_approximate(Φ::PeriodicKernelTranslates) =
    kernel_support_approximate(kernel(Φ))
kernel_support_approximate(kernel) = kernel_support(kernel)



## Application support

"Kernel that represents the derivative of another kernel."
struct DiffKernel{N,K}
    kernel  ::  K
    order   ::  NTuple{N,Int}
end

kernel_eval(k::DiffKernel{1}, x) = kernel_eval_derivative(k.kernel, x, k.order[1])
