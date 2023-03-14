module RefinableFunctions

using ..BasisTranslates
using BasisTranslates.BSplines

using BasisFunctions, DomainSets, CompositeTypes
using LinearAlgebra, BlockArrays

import DomainSets:
    numtype

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
    PeriodicRefinables,
    isrefinable,
    refinable_coeff,
    refinable_quadrature,
    cdfdual

# from util/sequences.jl
export datavector

include("util/common.jl")
include("util/conv.jl")
include("util/sequences.jl")
include("util/filterbank.jl")
include("inner.jl")
include("kernel.jl")
include("dyadic.jl")
include("dict.jl")
include("mra.jl")
include("wavelets.jl")
include("filters.jl")
include("scaling.jl")
include("splines.jl")

end
