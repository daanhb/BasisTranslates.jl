module RadialBasisFunctions

using ..BasisTranslates

using BasisFunctions, DomainSets
using LinearAlgebra

export MQ,
    IMQ,
    IQ,
    GA,
    PHS,
    TPS,
    PeriodicRBFs,
    PeriodicGaussians

import BasisTranslates:
    kernel,
    kernel_eval,
    kernel_eval_derivative,
    kernel_support,
    kernel_support_approximate

include("rbf_defs.jl")
include("kernel.jl")
include("dict.jl")

end
