module RationalFunctions

using ..BasisTranslates

using BasisFunctions, DomainSets, DomainIntegrals
using LinearAlgebra

import BasisTranslates:
    kernel,
    kernel_eval,
    kernel_eval_derivative,
    kernel_support,
    kernel_support_approximate

export PartialFractions,
    poles,
    lightning_uniform_poles,
    lightning_tapered_poles,
    lightning_samples

include("kernel.jl")
include("dict.jl")
include("lightning.jl")
include("projection.jl")
include("barycentric.jl")

end
