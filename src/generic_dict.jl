
"Dictionary of translates of a given kernel function."
struct KernelTranslates{S,T,K} <: SimpleTranslates{S,T}
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

centers(Φ::KernelTranslates) = Φ.centers

support(Φ::KernelTranslates) = Φ.support


"Dictionary of periodized translates of a given kernel function."
struct PeriodicKernelTranslates{S,T,K} <: UnitPeriodicTranslates{S,T}
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

parent_kernel(Φ::PeriodicKernelTranslates) = Φ.kernel

Base.size(Φ::PeriodicKernelTranslates) = (Φ.n,)

BasisTranslates.linearscaling(Φ::PeriodicKernelTranslates) = Φ.scaling
