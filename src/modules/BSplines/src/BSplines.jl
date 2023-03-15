module BSplines

using ..BasisTranslates

using BasisFunctions, DomainSets

using SpecialFunctions

import DomainSets:
    numtype

export eval_bspline,
    eval_bspline_derivative,
    eval_centered_bspline,
    eval_centered_bspline_derivative,
    eval_periodic_bspline,
    eval_periodic_bspline_derivative,
    eval_periodic_centered_bspline,
    eval_periodic_centered_bspline_derivative,
    RegularBSplines,
    PeriodicBSplines,
    spline_order,
    spline_degree,
    BSpline

include("bspline_defs.jl")
include("integral.jl")
include("kernel.jl")
include("dict.jl")

end
