module BasisTranslates

using BasisFunctions, DomainSets

# From utilities
export PeriodicInterval

# From the generic code
export translates_grid,
    translate_center

# utilities
include("util/periodicinterval.jl")

include("dictionary.jl")
include("periodic_dict.jl")

# submodules
include("modules/BSplines/BSplines.jl")
include("modules/RadialBasisFunctions/RadialBasisFunctions.jl")

end
