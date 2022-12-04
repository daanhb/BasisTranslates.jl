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

import BasisFunctions:
    coefficients

export RefinableFunction

include("util/conv.jl")
include("sequences.jl")
include("cache.jl")
include("kernel.jl")
include("dyadic.jl")
include("dict.jl")

end
