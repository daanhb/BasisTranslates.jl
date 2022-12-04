
"The Haar MRA."
struct HaarMRA{T} <: MRA{T}
end

HaarMRA() = HaarMRA{Float64}()

isorthogonal(::HaarMRA) = true

primal_lowpass(mra::HaarMRA{T}) where {T} =
    VectorSequence(1/sqrt(T(2))*[1,1])
