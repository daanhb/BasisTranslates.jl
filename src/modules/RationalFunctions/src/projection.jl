function rat_projection_matrix(basis, f, μ, ::Type{T}) where T
    N = length(basis)
    A = zeros(T,N,N)
    B = zeros(T,N,N)
    C = zeros(T,N,N)
    D = zeros(T,N,N)
    for i in 1:N
        for j in 1:N
            A[i,j] = integral(x->basis[i](x)*basis[j](x)*f(x)^2, μ)
            B[i,j] = -integral(x->2x*basis[i](x^2)*basis[j](x^2)*f(x^2), μ)
            C[i,j] = B[i,j]
            D[i,j] = integral(x->basis[i](x)*basis[j](x), μ)
        end
    end
    [A B; C D]
end

function smallest_singular_value_rational(basis, M)
    N = length(basis)
    u,s,v = svd(M)
    @show s
    U = u[:,end]
    V = v[:,end]
    p = Expansion(basis, V[1:N])
    q = Expansion(basis, V[N+1:2N])
    p, q
end

function rat_projection(f, N, domain::AbstractInterval{T}) where T
    basis = ChebyshevT{T}(N) → domain
    m = mapto(ChebyshevInterval{T}(), domain)
    μ = mappedmeasure(m, LegendreWeight{T}())
    M = rat_projection_matrix(basis, f, μ, T)
    p, q = smallest_singular_value_rational(basis, M)
    p, q
end
