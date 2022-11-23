## Integrals of B-splines

"""
Calculate ∫\$[B_n(x)]^2\$ dx, with \$B_n\$ the B-spline of degree n.
"""
squared_spline_integral(n::Int) = squared_spline_integral(n, BigInt)

function squared_spline_integral(N::Int, ::Type{T}) where {T<:AbstractFloat}
    m = convert(T,N+1)
    S = zero(T)
    for k in zero(T):m-1
        S += my_binomial(2m,k)*(-1)^k*(2*(m-k))^(2m-1)
    end
    S/(convert(T,2)^(2m-1)*my_factorial(2m-1))
end

function squared_spline_integral(N::Int, ::Type{T}) where {T<:Integer}
    m = convert(T, N+1)
    if m > 7 && sizeof(T) <= 8
        error("Higher precision is required, try `squared_spline_integral($N, BigInt)`")
    end
    S = zero(Rational{T})
    for k in 0:m-1
        S += my_binomial(2m,k)*(-1)^k*(2*(m-k))^(2m-1)
    end
    S//(2^(2m-1)*my_factorial(2m-1))
end


"""
Calculate ∫ \$B_m(x) B_m(x-k)\$ dx, with \$B_m\$ the B-spline of degree m.
"""
shifted_spline_integral(m::Int, k::Int) = shifted_spline_integral(m, k, BigInt)

function shifted_spline_integral(m::Int, k::Int, ::Type{T}) where {T<:AbstractFloat}
    @assert k >= 0
    E = convert(T, k)
    M = convert(T, m)
    I = zero(T)
    O = one(T)
    if k > m
        return zero(T)
    end
    for l in zero(T):M+O
        for j in zero(T):M+O
            S = zero(T)
            if j <= l-E
                for i in zero(T):M
                    S += my_binomial(M,i)*(l-j-E)^(m-i)*(M+O-l)^(i+M+O)/(i+m+O)
                end
                I += (-O)^(j+l)*my_binomial(M+O,l)*my_binomial(M+O,j)*S
            else
                for i in zero(T):M
                    S += my_binomial(M,i)*(E+j-l)^(M-i)*(M+O-j-E)^(i+M+O)/(i+M+O)
                end
                I += (-O)^(j+l)*my_binomial(M+O,l)*my_binomial(M+O,j)*S
            end
        end
    end
    convert(T,I/(my_factorial(M)^convert(T,2)))
end

function shifted_spline_integral(m::Int, k::Int, ::Type{T}) where {T<:Integer}
    k = abs(k)
    if k == 0
        return squared_spline_integral(m, T)
    end
    if k > m
        return convert(Rational{T},0)
    end
    m = convert(T,m)
    I = convert(Rational{T},0)
    if m > 7 && sizeof(T) <= 8
        error("You need a higher precision, try `shifted_spline_integral($m, $t, BigInt)`")
    end
    for l in 0:m+1
        for j in 0:m+1
            if j <= l-k
                S = zero(Rational{T})
                for i in 0:m
                    S += my_binomial(m,i)*(l-j-k)^(m-i)*(m+1-l)^(i+m+1)//(i+m+1)
                end
                I = I + convert(Rational{T},-1)^(j+l)*my_binomial(m+1,l)*my_binomial(m+1,j)*S
            else
                S = zero(Rational{T})
                for i in 0:m
                    S += my_binomial(m,i)*(k+j-l)^(m-i)*(m+1-j-k)^(i+m+1)//(i+m+1)
                end
                I = I + convert(Rational{T},-1)^(j+l)*my_binomial(m+1,l)*my_binomial(m+1,j)*S
            end
        end
    end
    I//(my_factorial(m)^2)
end
