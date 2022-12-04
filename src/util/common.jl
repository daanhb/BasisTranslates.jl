
const RegularGrid{T} = StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T},Int}

# This implementation is provided pending the resolution of #47717 in julia
# see: https://github.com/JuliaLang/julia/issues/47717
#
# function LinearAlgebra.pinv(D::Diagonal{T}; atol::Real = 0.0, rtol::Real = (eps(real(float(oneunit(T))))*length(D.diag))*iszero(atol)) where T
function diagonal_pinv(D::Diagonal{T}; atol::Real = 0.0, rtol::Real = (eps(real(float(oneunit(T))))*length(D.diag))*iszero(atol)) where T
    Di = similar(D.diag, typeof(inv(oneunit(T))))
    if !isempty(D.diag)
        maxabsA = maximum(abs, D.diag)
        abstol = max(rtol * maxabsA, atol)
        for i in 1:length(D.diag)
            if abs(D.diag[i]) > abstol
                invD = inv(D.diag[i])
                if isfinite(invD)
                    Di[i] = invD
                    continue
                end
            end
            # fallback
            Di[i] = zero(T)
        end
    end
    Diagonal(Di)
end
