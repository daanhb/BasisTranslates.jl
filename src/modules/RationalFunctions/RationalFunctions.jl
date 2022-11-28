module RationalFunctions

using ..BasisTranslates

using BasisFunctions, DomainSets

import BasisTranslates:
    kernel_eval,
    kernel_eval_derivative,
    kernel_support,
    kernel_support_approximate

export PartialFractions,
    poles,
    uniform_poles,
    tapered_poles

"The L1 norm."
struct L1Norm
end

include("partialfractions.jl")
include("lightning.jl")

end
