
"""
A dictionary of periodic radial basis functions on the interval `[0,1]`.
"""
abstract type PeriodicRBFs{T} <: BasisTranslates.PeriodicTranslates{T,T} end

Base.size(Φ::PeriodicRBFs) = (Φ.n,)

"What is the shape parameter of the kernels?"
shape_parameter(Φ::PeriodicRBFs) = Φ.epsilon

BasisTranslates.linearscaling(Φ::PeriodicRBFs) = false

function BasisTranslates.translate_map(Φ::PeriodicRBFs, idx)
    ε = shape_parameter(Φ)
    AffineMap(ε, -(idx-1)*ε/length(Φ))
end

# kernel_support_approximate(Φ::PeriodicRBFs, threshold = eps(prectype(Φ))) =
    # kernel_support_approximate(rbf_kernel(Φ), threshold)



"A dictionary of periodic Gaussian RBFs."
struct PeriodicGaussians{T} <: PeriodicRBFs{T}
    epsilon ::  T
    n       ::  Int
end

BasisTranslates.parent_kernel(Φ::PeriodicGaussians) = Gaussian()

"Optimal value of the linear scaling constant for Gaussian approximations."
optimal_ga_scaling(period, τ) = pi / (period*sqrt(2*log(1+τ^(-2))))

# kernel_support(Φ::PeriodicGaussians{T}) where {T} = FullSpace{T}()
