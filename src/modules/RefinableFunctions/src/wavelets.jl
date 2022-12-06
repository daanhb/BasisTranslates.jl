
"A basis of periodic wavelets."
struct PeriodicWavelets{T} <: BasisTranslates.Translates{T,T,:unitperiodic}
    wavelet     ::  CompactWavelet{T}
    n           ::  Int
end

Base.size(Φ::PeriodicWavelets) = (Φ.n,)
parent_kernel(Φ::PeriodicWavelets) = Φ.wavelet
Base.similar(Φ::PeriodicWavelets{T}, ::Type{T}, n::Int) where {T} =
    PeriodicWavelets(Φ.wavelet, n)



## A basis of wavelets

compute_offsets(mra, n, levels) = vcat([n>>levels], [n>>i for i in levels:-1:1])

function compute_dicts(mra::MRA{T}, n, levels) where {T}
    φ = EvaluatedRefinable(scalingfunction(mra))
    ψ = CompactWavelet(φ, primal_highpass(mra))
    n0 = n >> levels
    dict0 = PeriodicRefinables(φ, n0)
    wdicts = [PeriodicWavelets(ψ, n>>i) for i in levels:-1:1]
    [dict0; wdicts]
end

"A wavelet basis."
struct WaveletBasis{T} <: BasisFunctions.CompositeDict{T,T}
    mra         ::  MRA{T}
    n           ::  Int
    levels      ::  Int
    offsets     ::  Vector{Int}
    dicts       ::  Vector{Dictionary}

    function WaveletBasis{T}(mra::MRA, n::Int, levels::Int) where {T}
        offsets = compute_offsets(mra, n, levels)
        dicts = compute_dicts(mra, n, levels)
        new(mra, n, levels, offsets, dicts)
    end
end

WaveletBasis(mra::MRA{T}, args...) where {T} =
    WaveletBasis{T}(mra, args...)

Base.size(Φ::WaveletBasis) = (Φ.n,)
