module RefinableFunctions

using ..BasisTranslates

using BasisFunctions, DomainSets
using LinearAlgebra

import BasisTranslates:
    iscompact,
    support,
    kernel,
    kernel_eval,
    kernel_eval_derivative,
    kernel_support,
    kernel_support_approximate

export RefinableFunction

include("sequences.jl")
include("kernel.jl")
include("dyadic.jl")
include("dict.jl")

end
