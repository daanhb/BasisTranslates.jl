
"Supertype of discrete wavelets associated with a multi-resolution analysis."
abstract type DiscreteWavelet{T} <: BasisTranslates.Kernel end

"A compactly supported wavelet."
struct CompactWavelet{T} <: DiscreteWavelet{T}
    scalingfunction ::  Refinable{T}
    coefficients    ::  VectorSequence{T}
end

eval_kernel(ψ::CompactWavelet, x) = _eval_kernel(ψ, x, ψ.scalingfunction, ψ.coefficients)
_eval_kernel(ψ::CompactWavelet, x, φ, coef) =
    sum(coef[i]*eval_kernel(φ, 2x-i) for i in support(coef))


"Form the sequence `(-1)^k h_{n-k}` from the given sequence."
function alternating_flip(h::CompactSequence, n::Int = 1)
    @assert isodd(n)
    I = support(h)
    I2 = flipunitrange(I) .+ n
    VectorSequence([adjoint((-1)^(i+1)*h[i]) for i in reverse(I)], I2)
end


"""
Supertype of a multi-resolution analysis.

A multi-resolution analysis consists of a primal pair and a dual pair
of low-pass and high-pass filters.
"""
abstract type MRA{T} end

isorthogonal(mra::MRA) = false
isbiorthogonal(mra::MRA) = true

"Return the primal low-pass filter of the MRA."
function primal_lowpass(mra::MRA) end

"Return the dual low-pass filter of the MRA."
dual_lowpass(mra::MRA) =
    isorthogonal(mra) ? primal_lowpass(mra) : error("Don't know dual low-pass filter of MRA.")

"Return the primal high-pass filter of the MRA."
primal_highpass(mra::MRA) = alternating_flip(dual_lowpass(mra))

"Return the dual high-pass filter of the MRA."
dual_highpass(mra::MRA) = alternating_flip(primal_lowpass(mra))

function analysis_lowpass(mra::MRA, n::Int)
    lo_d = dual_lowpass(mra)
    periodic_downsampled_filter_array(adjoint(lo_d), n)
end
function analysis_highpass(mra::MRA, n::Int)
    hi_d = dual_highpass(mra)
    periodic_downsampled_filter_array(adjoint(hi_d), n)
end
function synthesis_lowpass(mra::MRA, n::Int)
    lo_r = primal_lowpass(mra)
    periodic_upsampling_filter_array(lo_r, n)
end
function synthesis_highpass(mra::MRA, n::Int)
    hi_r = primal_highpass(mra)
    periodic_upsampling_filter_array(hi_r, n)
end

const AnalysisMRA{T} =
    BlockMatrix{T, Matrix{FatRectangularCirculant{T}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

const SynthesisMRA{T} =
    BlockMatrix{T, Matrix{TallRectangularCirculant{T}}, Tuple{BlockedUnitRange{Vector{Int64}}, BlockedUnitRange{Vector{Int64}}}}

function periodic_analysis(mra::MRA, n::Int)
    A1 = analysis_lowpass(mra, n)
    A2 = analysis_highpass(mra, n)
    mortar(reshape([A1,A2],2,1))
end

function periodic_synthesis(mra::MRA, n::Int)
    B1 = synthesis_lowpass(mra, n)
    B2 = synthesis_highpass(mra, n)
    mortar(reshape([B1,B2],1,2))
end
