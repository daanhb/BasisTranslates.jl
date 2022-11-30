module RationalFunctions

using ..BasisTranslates

using BasisFunctions, DomainSets

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

include("partialfractions.jl")
include("lightning.jl")

end
