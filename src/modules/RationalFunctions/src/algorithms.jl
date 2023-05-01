"""
Compute the matrix corresponding to the linearized projection method.

This matrix is comparable to a Gram matrix of a continuous least squares problem.
"""
function rat_projection_matrix(basis, f, μ, ::Type{T}) where T
    N = length(basis)
    A = zeros(T,N,N)
    B = zeros(T,N,N)
    C = zeros(T,N,N)
    D = zeros(T,N,N)
    for i in 1:N
        for j in 1:N
            A[i,j] = integral(x->basis[i](x)*basis[j](x)*f(x)^2, μ)
            B[i,j] = -integral(x->basis[i](x)*basis[j](x)*f(x), μ)
            C[i,j] = B[i,j]
            D[i,j] = integral(x->basis[i](x)*basis[j](x), μ)
        end
    end
    [A B; C D]
end

function smallest_singular_value_rational(basis, M)
    N = length(basis)
    u,s,v = svd(M)
    V = v[:,end]
    p = Expansion(basis, V[1:N])
    q = Expansion(basis, V[N+1:2N])
    p, q
end

"""
Compute a degree `N` rational approximation by linearized projection.

The method computes the polynomials `p` and `q` that minimize `||pf-q||`,
which is a continuous norm on the given domain.
"""
function rat_projection(f, N, domain::AbstractInterval{T}) where T
    m = mapto(ChebyshevInterval{T}(), domain)
    μ = mappedmeasure(m, LegendreWeight{T}())
    rat_projection(f, N, μ)
end

function rat_projection(f, N, μ::Measure{T}) where T
    basis = ChebyshevT{T}(N) → support(μ)
    M = rat_projection_matrix(basis, f, μ, T)
    p, q = smallest_singular_value_rational(basis, M)
    p, q, M
end


"""
Compute the AAA rational approximation to the function values `F`
in the points `X`.
"""
function aaa(F, X; tol = 1e-13, degree = -1)
    local wj, α, β

    XT = eltype(X)
    FT = eltype(F)
    Z = X
    M = length(Z)
    @assert length(F) == M

    reltol = tol * norm(F, Inf)
    SF = Diagonal(F)
    R = sum(F)/M*ones(FT, size(F))
    mmax = 100
    if degree > 0
        mmax = degree+1
    end

    zj = XT[]
    fj = FT[]
    C = zeros(FT,M,0)
    J = 1:M
    errvec = real(XT)[]

    for m in 1:mmax
        # Select next support point where error is largest:
        err,jj = findmax(abs.(F-R))
        push!(zj, Z[jj])
        push!(fj, F[jj])
        J = setdiff(J, jj)
        C = [C 1 ./ (Z .- Z[jj])]

        # Compute weights
        Sf = Diagonal(fj)
        A = SF*C - C*Sf
        if length(J) > 0
            U,S,V = svd(A[J,:])
            wj = V[:,m]
        else
            wj = zeros(XT,m)
            wj[end] = 1
        end

        α = wj .* fj
        β = wj

        # Rational approximant on Z:
        N = C*α
        D = C*β
        R = copy(F)
        R[J] = N[J] ./ D[J]

        # Error in the sample points:
        maxerr = norm(F-R, Inf)
        errvec = [errvec; maxerr]

        # Check if converged
        if maxerr < reltol
            break
        end
    end

    r = z -> barycentric_eval(zj, α, β, z)
    pol, res, zer = barycentric_poles_residues_zeros(zj, α, β)
    r, pol, res, zer, zj, fj, wj, errvec
end

"Apply the AAA algorithm and return the result in barycentric form."
function barycentric_aaa(F, X; options...)
    r, pol, res, zer, zj, fj, wj, errvec = aaa(F, X; options...)
    α = wj .* fj
    β = wj
    Expansion(BarycentricRational(zj, β), α)
end
