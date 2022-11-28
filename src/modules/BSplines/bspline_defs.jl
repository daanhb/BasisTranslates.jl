# The code in this file is based on the initial implementation
# by vincentcp in an unregistered package called CardinalBSplines.jl.

_prectype(x) = promote_type(Float64, typeof(x))
_prectype(x,y) = promote_type(Float64, typeof(x), typeof(y))

my_binomial(n, k) = Base.binomial(n, k)
my_binomial(n::T, k::T) where {T <: AbstractFloat} = gamma(n+1)/gamma(k+1)/gamma(n-k+1)

my_factorial(n) = Base.factorial(n)
my_factorial(n::T) where {T<:AbstractFloat} = gamma(n+1)

function periodize(x, period)
    x -= period*fld(x, period)
    x >= period && (x -= period)
    # @assert(0<= x < period)
    x
end


support_bspline(n::Int, ::Type{T}) where {T} = (zero(T), T(n+1))
support_centered_bspline(n::Int, ::Type{T}) where {T} = (-T(n+1)/2, T(n+1)/2)



"""
Evaluate the periodic centered B-spline of degree `n`, with the
given `period`, in the point `x`.
"""
eval_periodic_centered_bspline(n::Int, x, period, ::Type{T} = _prectype(x,period)) where {T} =
    eval_periodic_bspline(n, x+T(n+1)/2, period, T)

"""
Evaluate the centered B-spline of degree `n` in the point `x`.

The support of the centered B-spline is \$\\left[-\\tfrac{n+1}{2},\\tfrac{n+1}{2}\\right]\$.
"""
eval_centered_bspline(n::Int, x, ::Type{T} = _prectype(x)) where {T} =
    eval_bspline(n, x+T(n+1)/2, T)

"""
Evaluate the periodic B-spline of degree `n`, with the given `period`,
in the point 'x'.
"""
function eval_periodic_bspline(n::Int, x, period, ::Type{T} = _prectype(x,period)) where {T}
    x = periodize(x, period)
    sum(eval_bspline(n, x+period*k, T) for k in 0:floor(Int, (n+1-x)/period); init = zero(T))
end

"""
Evaluate the B-spline of degree `n` in the point `x`.

The support of the B-spline is `[0,n+1]`.
"""
function eval_bspline(n::Int, x, ::Type{T} = _prectype(x)) where {T}
    if n == 0
        eval_bspline_degree0(x, T)
    elseif n == 1
        eval_bspline_degree1(x, T)
    elseif n == 2
        eval_bspline_degree2(x, T)
    elseif n == 3
        eval_bspline_degree3(x, T)
    elseif n == 4
        eval_bspline_degree4(x, T)
    else
        @assert n > 4
        (T(x)/T(n)*eval_bspline(n-1, x, T) +
            (T(n+1)-T(x))/T(n)*eval_bspline(n-1, x-1, T))
    end
end

eval_bspline_degree0(x, ::Type{T}) where {T} = (0 <= x < 1) ? one(T) : zero(T)

function eval_bspline_degree1(x, ::Type{T}) where {T}
    if (0 <= x < 1)
        return T(x)
    elseif (1 <= x < 2)
        return T(2) - T(x)
    else
        return T(0)
    end
end

function eval_bspline_degree2(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x),T(0), T(0), T(1/2))
    elseif (1 <= x < 2)
        return @evalpoly(T(x),T(-3/2), T(3), T(-1))
    elseif (2 <= x < 3)
        return @evalpoly(T(x),T(9//2), T(-3), T(1//2))
    else
        return T(0)
    end
end

function eval_bspline_degree3(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(0), T(1//6))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(2//3), T(-2), T(2), T(-1//2))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-22//3), T(10), T(-4), T(1//2))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(32//3), T(-8), T(2), T(-1//6))
    else
        return T(0)
    end
end

function eval_bspline_degree4(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(0), T(0), T(1//24))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(-5//24), T(5//6), T(-5/4), T(5//6), T(-1//6))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(155//24), T(-25//2), T(35//4), T(-5//2), T(1//4))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(-655//24), T(65//2), T(-55//4), T(5//2), T(-1//6))
    elseif (4 <= x < 5)
        return @evalpoly(T(x), T(625//24), -T(125//6), T(25//4), T(-5//6), T(1//24))
    else
        return T(0)
    end
end

"""
Evaluate the `d`th derivative of the B-spline of degree `n`,
in the point `x`.
"""
function eval_bspline_derivative(n::Int, x, order::Int, ::Type{T} = _prectype(x)) where {T}
    @assert n >= 0
    @assert order >= 0
    if order > n
        zero(T)
    elseif order == 0
        eval_bspline(n, x, T)
    elseif order == n
        xfl = floor(Int, x)
        if xfl < 0 || xfl > n
            T(0)
        else
            T((-1)^xfl*my_binomial(n,xfl))
        end
    else
        # reduce degree
        T(d)/T(n)*(eval_bspline_derivative(n-1, x, order-1, T) - eval_bspline_derivative(n-1, x-1, order-1, T)) +
            (T(n+1)-T(x))/T(n)*eval_bspline_derivative(n-1, x-1, order, T) + T(x)/T(n)*eval_bspline_derivative(n-1, x, order, T)
    end
end

eval_bspline_derivative_degree0(d::Int, x, ::Type{T} = _prectype(x)) where {T} =
    d == 0 ? eval_bspline(0, x, T) : zero(T)

function eval_bspline_derivative_degree1_diff1(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return - T(1)
    else
        return T(0)
    end
end

function eval_bspline_derivative_degree2_diff1(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x),T(0), T(1))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(3), T(-2))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-3), T(1))
    else
        return T(0)
    end
end

function eval_bspline_derivative_degree3_diff1(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(1//2))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(-2), T(4), T(-3//2))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(10), T(-8), T(3//2))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(-8), T(4), T(-1//2))
    else
        return T(0)
    end
end

function eval_bspline_derivative_degree4_diff1(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(0), T(1//6))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(5//6), T(-5/2), T(5//2), T(-2//3))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-25//2), T(35//2), T(-15//2), T(1))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(65//2), T(-55//2), T(15//2), T(-2//3))
    elseif (4 <= x < 5)
        return @evalpoly(T(x), -T(125//6), T(25//2), T(-5//2), T(1//6))
    else
        return T(0)
    end
end


function eval_bspline_derivative_degree2_diff2(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return T(-2)
    elseif (2 <= x < 3)
        return T(1)
    else
        return T(0)
    end
end

function eval_bspline_derivative_degree3_diff2(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(1))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(4), T(-3))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-8), T(3))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(4), T(-1))
    else
        return T(0)
    end
end

function eval_bspline_derivative_degree4_diff2(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(1//2))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(-5/2), T(5), T(-2))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(35//2), T(-15), T(3))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(-55//2), T(15), T(-2))
    elseif (4 <= x < 5)
        return @evalpoly(T(x), T(25//2), T(-5), T(1//2))
    else
        return T(0)
    end
end

function eval_bspline_derivative_degree3_diff3(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return T(-3)
    elseif (2 <= x < 3)
        return  T(3)
    elseif (3 <= x < 4)
        return T(-1)
    else
        return T(0)
    end
end

function eval_bspline_derivative_degree4_diff3(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(1))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(5), T(-4))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-15), T(6))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(15), T(-4))
    elseif (4 <= x < 5)
        return @evalpoly(T(x), T(-5), T(1))
    else
        return T(0)
    end
end

function eval_bspline_derivative_degree4_diff4(x, ::Type{T} = _prectype(x)) where {T}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return T(-4)
    elseif (2 <= x < 3)
        return T(6)
    elseif (3 <= x < 4)
        return T(-4)
    elseif (4 <= x < 5)
        return  T(1)
    else
        return T(0)
    end
end


"""
Evaluate the `d`th derivative of the periodic B-spline of degree `n` in `x`.
"""
function eval_periodic_bspline_derivative(n::Int, x, period, order::Int, ::Type{T} = _prectype(x,period)) where {T}
    x = periodize(x, period)
    sum(eval_bspline_derivative(n, x+period*k, order, T) for k in 0:floor(Int, (n+1-x)/period); init=zero(T))
end

"""
Evaluate the `d`th derivative of the periodic centered B-spline of degree `n`,
with the given `period`, in the point `x`.
"""
eval_periodic_centered_bspline_derivative(n::Int, x, period, order::Int, ::Type{T} = _prectype(x,period)) where {T} =
    eval_periodic_bspline_derivative(n, x+T(n+1)/2, period, order, T)


"""
Evaluate the `d`th derivative of the centered B-spline of degree `n`,
with the given `period`, in the point `x`.
"""
eval_centered_bspline_derivative(n::Int, x, order::Int, ::Type{T} = _prectype(x)) where {T} =
    eval_bspline_derivative(n, x+T(n+1)/2, order, T)
