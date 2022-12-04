module RefinableFunctions

using ..BasisTranslates
using BasisTranslates.BSplines

using BasisFunctions, DomainSets
using LinearAlgebra, BlockArrays

import BasisTranslates:
    iscompact,
    support,
    offset,
    kernel,
    parent_kernel,
    kernel_eval,
    kernel_eval_derivative,
    kernel_support,
    coefficients,
    isorthogonal,
    isbiorthogonal

export RefinableFunction,
    PeriodicRefinables

include("util/common.jl")
include("util/conv.jl")
include("util/sequences.jl")
include("util/filterbank.jl")
include("kernel.jl")
include("dyadic.jl")
include("dict.jl")
include("mra.jl")
include("filters.jl")
include("scaling.jl")

end
