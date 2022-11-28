
"""
A dictionary of periodic radial basis functions on the interval `[0,1]`.
"""
abstract type PeriodicRBFs{T} <: BasisTranslates.PeriodicTranslates{T,T} end

Base.size(Φ::PeriodicRBFs) = (Φ.n,)

"What is the shape parameter of the kernels?"
shape_parameter(Φ::PeriodicRBFs) = Φ.epsilon

kernel_support_approximate(Φ::PeriodicRBFs, threshold = eps(prectype(Φ))) =
    support_approximate(rbf_kernel(Φ), threshold) / shape_parameter(Φ)

kernel_eval(Φ::PeriodicRBFs, x) =
    _kernel_eval(rbf_kernel(Φ), shape_parameter(Φ), x)
_kernel_eval(rbf::RBF, epsilon, x) = kernel_eval(rbf, epsilon*x)

kernel_eval_derivative(Φ::PeriodicRBFs, order, x) =
    _kernel_eval_derivative(rbf_kernel(Φ), shape_parameter(Φ), order, x)
_kernel_eval_derivative(rbf::RBF, epsilon, order, x) =
    epsilon^order*kernel_eval_derivative(rbf, order, epsilon*x)

BasisTranslates.linearscaling(Φ::PeriodicRBFs) = false


"A dictionary of periodic Gaussian RBFs."
struct PeriodicGaussians{T} <: PeriodicRBFs{T}
    epsilon ::  T
    n       ::  Int
end

rbf_kernel(Φ::PeriodicGaussians) = Gaussian()

"Optimal value of the linear scaling constant for Gaussian approximations."
optimal_ga_scaling(n, τ) = 2pi / sqrt(2*pi*log(1+τ^(-2)))
