
function db2_filter(::Type{T} = Float64) where {T}
    sq2 = sqrt(T(2))
    sq3 = sqrt(T(3))
    coef = 1/sq2*[(1+sq3)/4, (3+sq3)/4, (3-sq3)/4, (1-sq3)/4]
    VectorSequence(coef)
end
