module BasisTranslates

using BlockArrays, LinearAlgebra
using BasisFunctions, DomainSets
using AbstractFFTs

# From utilities
export PeriodicInterval

# From the generic code
export translates_grid,
    translate_center

# utilities
include("util/common.jl")
include("util/periodicinterval.jl")
include("util/fourier.jl")
include("util/permutations.jl")
include("util/blockarrays.jl")
include("util/multicirculant.jl")

include("dictionary.jl")
include("periodic_dict.jl")

# submodules
include("modules/BSplines/BSplines.jl")
include("modules/RadialBasisFunctions/RadialBasisFunctions.jl")

end
