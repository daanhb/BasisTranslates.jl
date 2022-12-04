
"The Haar MRA."
struct HaarMRA{T} <: MRA{T}
end

HaarMRA() = HaarMRA{Float64}()

isorthogonal(::HaarMRA) = true

primal_lowpass(mra::HaarMRA{T}) where {T} =
    VectorSequence(1/sqrt(T(2))*[1,1])


"The Daubechies family MRA."
struct DaubechiesMRA{T} <: MRA{T}
    order   ::  Int
end

DaubechiesMRA(order::Int) = DaubechiesMRA{Float64}(order)

isorthogonal(::DaubechiesMRA) = true

function primal_lowpass(mra::DaubechiesMRA{T}) where {T}
    p = mra.order
    if p == 1
        primal_lowpass(HaarMRA{T}())
    elseif p == 2
        db2(T)
    else
        error("Higher order Daubechies scaling functions not yet implemented.")
    end
end
