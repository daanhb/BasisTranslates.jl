
uniform_poles(σ::Int, n) = uniform_poles(float(σ), n)
uniform_poles(σ::T, n) where {T} = -exp.(-σ*(0:n-1)./sqrt(T(n)))

tapered_poles(σ::Int, n) = tapered_poles(float(σ), n)
tapered_poles(σ::T, n) where {T} = -exp.(σ*(sqrt.(T(n):-1:1) .- sqrt(T(n))))

lightning_normalization(Φ::PartialFractions) =
    DiagonalOperator(Φ, -Φ.poles)
