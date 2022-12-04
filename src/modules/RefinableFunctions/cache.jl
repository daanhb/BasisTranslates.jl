
"Supertype of caches to store evaluation of functions."
abstract type EvaluationCache end

const PointRange{T} = StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T},Int}

"Store function evaluations in equispaced points."
struct EquispacedCache{T} <: EvaluationCache
    pts     ::  PointRange{T}
    vals    ::  Vector{T}
end

# Perform linear interpolation on the function samples
function eval(cache::EquispacedCache{T}, x) where {T}
    a = first(cache.pts)
    b = last(cache.pts)
    if a <= x <= b
        h = step(cache.pts)
        i = floor(Int, (x-a)/h)
        if i == length(cache.pts)
            cache.vals[end]
        else
            t1 = cache.pts[i]
            t2 = cache.pts[i+1]
            z1 = cache.vals[i]
            z2 = cache.vals[i+1]
            z1 + (x-t1)/(t2-t1)*(z2-z1)
        end
    else
        zero(T)
    end
end
