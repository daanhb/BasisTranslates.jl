
# Kernels do not need to inherit from Kernel, but it is convenient
# to provide a fallback interface.

"Supertype of kernels."
abstract type Kernel end

issymetric(φ::Kernel) = false

kernel_support_approximate(φ::Kernel, threshold...) = kernel_support(φ)

hascompactsupport(φ::Kernel) = iscompact(kernel_support(φ))
hascompactsupport_approximate(φ::Kernel, threshold...) =
    iscompact(kernel_support_approximate(φ, threshold...))


"Kernel that represents the derivative of another kernel."
struct DiffKernel{N,K} <: BasisTranslates.Kernel
    kernel  ::  K
    order   ::  NTuple{N,Int}
end

DiffKernel(φ, order::Int) = DiffKernel(φ, (order,))

parent_kernel(φ::DiffKernel) = φ.kernel

kernel_eval(φ::DiffKernel{1}, x) = kernel_eval_derivative(φ.kernel, x, φ.order[1])
kernel_support(φ::DiffKernel) = kernel_support(parent_kernel(φ))
kernel_support_approximate(φ::DiffKernel, threshold...) =
    kernel_support_approximate(parent_kernel(φ), threshold...)


"The kernel `φ(A*x)` for some factor `A`."
struct ScaledKernel{K,T} <: BasisTranslates.Kernel
    kernel  ::  K
    factor  ::  T
end

parent_kernel(φ::ScaledKernel) = φ.kernel

kernel_eval(φ::ScaledKernel, x) = kernel_eval(φ.kernel, φ.factor*x)
kernel_support(φ::ScaledKernel) = kernel_support(parent_kernel(φ))/φ.factor
kernel_support_approximate(φ::ScaledKernel, threshold...) =
    kernel_support_approximate(parent_kernel(φ), threshold...)/φ.factor
