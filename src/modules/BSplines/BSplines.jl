module BSplines

using ..BasisTranslates

using BasisFunctions, DomainSets

using SpecialFunctions

export eval_bspline,
    eval_bspline_derivative,
    eval_centered_bspline,
    eval_centered_bspline_derivative,
    eval_periodic_bspline,
    eval_periodic_bspline_derivative,
    eval_periodic_centered_bspline,
    eval_periodic_centered_bspline_derivative,
    PeriodicBSplines,
    splinedegree

include("bspline_defs.jl")
include("integral.jl")
include("kernel.jl")
include("dict.jl")

end
