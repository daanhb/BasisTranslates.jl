
lightning_uniform_poles(σ::Number, n) = lightning_uniform_poles(float(σ), n)
lightning_uniform_poles(σ::T, n) where {T <: AbstractFloat} =
    -exp.(-σ*(0:n-1)./sqrt(T(n)))

lightning_tapered_poles(σ::Number, n) = lightning_tapered_poles(float(σ), n)
lightning_tapered_poles(σ::T, n) where {T <: AbstractFloat} =
    -exp.(σ*(sqrt.(T(n):-1:1) .- sqrt(T(n))))

function lightning_samples(poles::AbstractVector{T}, osf::Int) where {T}
    n = length(poles)
    m = osf*n
    samples = zeros(T, m)
    pts = [0; sort(-poles)]
    for i in 1:n
        samples[(i-1)*osf+1:i*osf] = range(pts[i], stop=pts[i+1], length=osf+2)[2:end-1]
    end
    samples
end

lightning_normalization(Φ::BasisTranslates.Translates) =
    DiagonalOperator(Φ, -translates_grid(Φ))
