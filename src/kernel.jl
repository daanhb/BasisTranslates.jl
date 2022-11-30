
# Kernels do not need to inherit from Kernel, but it is convenient
# to provide a fallback interface.

"Supertype of kernels."
abstract type Kernel end

issymetric(φ::Kernel) = false

kernel_support_approximate(φ::Kernel, threshold...) = kernel_support(φ)

hascompactsupport(φ::Kernel) = iscompact(kernel_support(φ))
hascompactsupport_approximate(φ::Kernel, threshold...) =
    iscompact(kernel_support_approximate(φ, threshold...))
