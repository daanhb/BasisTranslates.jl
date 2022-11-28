
"Dictionary of translates of a given kernel function."
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
