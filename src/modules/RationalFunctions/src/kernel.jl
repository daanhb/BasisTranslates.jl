
"The `1/x` kernel."
struct PartialFractionKernel <: BasisTranslates.Kernel
end

kernel_eval(φ::PartialFractionKernel, x) = 1/x
kernel_support(φ::PartialFractionKernel) = FullSpace{Float64}()
