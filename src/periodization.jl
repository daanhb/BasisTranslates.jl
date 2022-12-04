
"Supertype of periodic kernels."
abstract type PeriodicKernel <: Kernel end

isperiodic(φ::Kernel) = false
isperiodic(φ::PeriodicKernel) = true

"Periodization of a parent kernel."
struct PeriodizedKernel{K,T} <: PeriodicKernel
    parent  ::  K
    period  ::  T
end

parent_kernel(φ::PeriodizedKernel) = φ.parent
period(φ::PeriodizedKernel) = φ.period

"""
    k_left, k_right = nb_overlapping_kernels(φ::Kernel, period)

The number of shifted kernels that overlap, due to periodicity. The functions
`φ(x + k*P)` may overlap with `φ(x)` for all `k` in `1:k_left`, and similarly
for `φ(x - k*P)` and `1:k_right`.
"""
function nb_overlapping_kernels(φ, period)
    a,b = extrema(kernel_support_approximate(φ))
    k = floor(Int, (b-a)/period)
    k, k
end

"The number of overlapping kernels at a point `x`."
function nb_overlapping_kernels(φ, period, x)
    a,b = extrema(kernel_support_approximate(φ))
    k_left = floor(Int, (b-x)/period)
    k_right = floor(Int, (x-a)/period)
    k_left, k_right
end

"Evaluate the periodization of the given kernel with the given period."
function periodized_kernel_eval(φ, period, x)
    k_left, k_right = nb_overlapping_kernels(φ, period, x)
    z = kernel_eval(φ, x)
    for k in 1:k_right
        z1 = kernel_eval(φ, x-k*period)
        z += z1
    end
    for k in 1:k_left
        z1 = kernel_eval(φ, x+k*period)
        z += z1
    end
    z
end

"Evaluate the derivative of the periodization of the given kernel."
function periodized_kernel_eval_derivative(φ, period, x, order)
    k_left, k_right = nb_overlapping_kernels(φ, period, x)
    z = kernel_eval_derivative(φ, x, order)
    for k in 1:k_right
        z1 = kernel_eval_derivative(φ, x-k*period, order)
        z += z1
    end
    for k in 1:k_left
        z1 = kernel_eval_derivative(φ, x+k*period, order)
        z += z1
    end
    z
end

kernel_eval(φ::PeriodizedKernel, x) =
    periodized_kernel_eval(parent_kernel(φ), period(φ), x)

kernel_eval_derivative(φ::PeriodizedKernel, x, order) =
    periodized_kernel_eval_derivative(parent_kernel(φ), period(φ), x, order)

kernel_support(φ::PeriodizedKernel) = kernel_support(parent_kernel(φ))

kernel_support_approximate(φ::PeriodizedKernel, args...) =
    kernel_support_approximate(parent_kernel(φ), args...)
