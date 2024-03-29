module BasisTranslates

using AbstractFFTs, FFTW
using BlockArrays, LinearAlgebra
using BasisFunctions, DomainSets, DomainIntegrals
using AZalgorithm

import BasisFunctions:
    isperiodic, period,
    support

# From utilities
export PeriodicInterval

# From the generic code
export centers,
    center,
    kernel,
    kernel_eval,
    kernel_eval_derivative,
    kernel_support,
    az_approximate,
    KernelTranslates,
    PeriodicKernelTranslates

# utilities
include("util/common.jl")
include("util/periodicinterval.jl")
include("util/domains.jl")
include("util/fourier.jl")
include("util/permutations.jl")
include("util/blockarrays.jl")
include("util/multicirculant.jl")

include("kernel.jl")
include("periodization.jl")
include("dictionary.jl")
include("periodic_dict.jl")
include("generic_dict.jl")
include("shifted.jl")

# submodules
include("modules/BSplines/src/BSplines.jl")
include("modules/RadialBasisFunctions/src/RadialBasisFunctions.jl")
include("modules/RationalFunctions/src/RationalFunctions.jl")
include("modules/RefinableFunctions/src/RefinableFunctions.jl")

end
