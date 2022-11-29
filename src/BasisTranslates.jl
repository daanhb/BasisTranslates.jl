module BasisTranslates

using AbstractFFTs, BlockArrays, LinearAlgebra, LinearMaps
using BasisFunctions, DomainSets

LinearMap = LinearMaps.LinearMap

# From utilities
export PeriodicInterval

# From the generic code
export translates_grid,
    translate_center,
    az_approximate,
    KernelTranslates,
    PeriodicKernelTranslates

# utilities
include("util/common.jl")
include("util/periodicinterval.jl")
include("util/fourier.jl")
include("util/permutations.jl")
include("util/blockarrays.jl")
include("util/multicirculant.jl")

include("az.jl")
include("dictionary.jl")
include("periodic_dict.jl")
include("kernel.jl")

# submodules
include("modules/BSplines/BSplines.jl")
include("modules/RadialBasisFunctions/RadialBasisFunctions.jl")
include("modules/RationalFunctions/RationalFunctions.jl")

end
