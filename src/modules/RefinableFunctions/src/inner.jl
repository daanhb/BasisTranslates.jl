
"A matrix with the convolution of two sequences on its rows, shifted by two on consecutive rows."
function transfer_matrix(s1::CompactSequence, s2::CompactSequence)
    # TODO: I am not sure this routine is correct
    s = convolve(s1, s2)
    i1 = first(support(s))
    i2 = last(support(s))
    d = datavector(s)
    L = length(d)
    L1 = length(support(s1))
    L2 = length(support(s2))
    L0 = max(L1,L2)
    A = zeros(L-2,L-2)
    for i in 1:L-2
        A[i,:] = s[i2-2i+1:i2-2i+L-2]
    end
    A
end

"Return the integrals `φ(x)*φ(x-k)` for all k."
function productintegral_moments(s::CompactSequence; squarednorm = 1)
    A = transfer_matrix(s, s)
    vals = eigenvector_with_eigenvalue_1(A)
    mom = VectorSequence(vals, -datalength(s)+2:datalength(s)-2)
    squarednorm/mom[0] * mom
end

function productintegral_moments(s1::CompactSequence, s2::CompactSequence; mixednorm = 1)
    A = transfer_matrix(s1, s2)
    vals = eigenvector_with_eigenvalue_1(A)
    # TODO: figure out where the range of the moments is
    VectorSequence(vals)
end

"""
Compute a discrete dual for the given kernel function using an oversampling
factor `ofs` and allowing `Q+1` coefficients.

This routine is experimental and possibly not entirely correct in all cases.
Please check the validity of the output, e.g., whether it respects symmetry
if applicable.
"""
function discrete_dual(kern, ofs, Q)
    S = kernel_support(kern)
    s1 = leftendpoint(S)
    s2 = rightendpoint(S)
    L = Int(s2-s1)
    A = [kernel_eval(kern, k/ofs-l) for l in -1:2L+1, k in 0:Q-1]
    B = zeros(size(A,1))
    B[L+2] = 1
    A\B
end

"""
Compute an oversampled quadrature rule associated with a primal and
dual function with a difference in levels.
"""
function oversampled_quadrature(kernel, dualkernel, level)
    osf = 2^level
    # first compute the quadrature points
    htilde = refinable_coeff(dualkernel)
    Idual = support(htilde)
    t1 = first(Idual)
    t2 = last(Idual)
    t = range(t1, stop=t2, step = 1/osf)
    m = length(t)
    # compute the basis functions that overlap with these points on the given level
    I = support(refinable_coeff(kernel))
    Q = 10
    skernel = BasisTranslates.ScaledKernel(kernel, 2.0^level)
    Phi = BasisTranslates.KernelTranslates(skernel, -Q:2.0^(-level):Q, -Q..Q)
    A = [2^(level/2)*Phi[k](t[i]) for k in eachindex(Phi), i in eachindex(t)]
    A_I = -Q*2^level:Q*2^level
    c = DiracSequence()
    for i in 1:level
        c = sum(htilde[i]*shift(c,i) for i in support(htilde))
    end
    B = [c[k] for k in A_I]
    sumA = sum(A; dims=2)[:]
    I1 = findfirst(sumA .> 0)
    I2 = findlast(sumA .> 0)
    A = A[I1:I2,:]
    B = B[I1:I2]
    w = A\B
    t, w
end


"""
Compute a quadrature rule for integrals involving a refinable function.

The method is described in:
W. Sweldens and R. Piessens, "Quadrature formulae and asymptotic error expansions
for wavelet approximations of smooth functions", SIAM J. Numer. Anal. 31(4),
pp. 1240-1264, 1994.

We return a quadrature rule with a fixed number `osf` of points per unit length.
This value is like an oversampling factor, hence its name.
"""
function refinable_quadrature(kernel::Kernel, osf = 1)
    coefficients = datavector(refinable_coeff(kernel))
    x, w = refinable_quadrature(coefficients, osf*length(coefficients)-1)
    w = w / sum(w) * refinable_moment(kernel)
    a,b = extrema(kernel_support(kernel))
    x = a .+ x * (b-a) / (length(coefficients)-1)
    x, w
end

function refinable_quadrature(coefficients::AbstractVector, M = length(coefficients)-1)
    @assert sum(coefficients) ≈ sqrt(2)
    L = length(coefficients)
    # Perform calculations in BigFloat for greater accuracy
    if isodd(L) || (iseven(L) && iseven(M))
        t = range(-BigFloat(1.),stop = BigFloat(1.), length = M+1)
    else
        # We have to shift the grid in order to include the center point
        t = range(-BigFloat(1.),stop = BigFloat(1.), length = 2M+1)[2:2:end-1]
    end
    A = chebyshev_vandermonde_matrix(t)
    m = chebyshev_moments(BigFloat.(coefficients), length(t)-1)
    w = A\m
    T = eltype(coefficients)
    map(T, (L-1) * (t .+ 1)/2), map(T, w)[:]
end

"Scale the given quadrature rule from `[a,b]` to `[c,d]`."
function scale_quadrature(x, w, a, b, c, d)
    L = b-a
    x2 = c .+ (x .- a)/L*(d-c)
    w2 = w * (d-c)/L
    x2, w2
end

function chebyshev_vandermonde_matrix(t::AbstractVector{T}) where {T}
    M = length(t)-1
    A = zeros(T, M+1, M+1)
    for i in 0:M
        A[i+1,:] = chebyshev_eval.(i, t)
    end
    A
end

"Evaluate the chebyshev polynomial of degree p in x."
chebyshev_eval(p, x) = real(cos(p*acos(x)))


"""
Calculate the modified moments of the scaling function defined by coefficients.
"""
function chebyshev_moments(coefficients::AbstractVector, M::Int)
    T = eltype(coefficients)
    m = zeros(T,M+1)
    m[1] = 1
    L = length(coefficients)

    w = zeros(T,M+1,M+1,L)
    for k in 1:L
        w[:,:,k] = wcoefficients(T, L, k, M)
    end

    for p in 1:M
        d = zero(T)
        for i in 0:p-1
            e = zero(T)
            for k in 1:L
                # % CAVE: our coefficients differ by a factor 1/sqrt(2) from Sweldens!
                e = e + 1/sqrt(T(2))*coefficients[k]*w[p+1,i+1,k]
            end
            d = d + e*m[i+1]
        end
        m[p+1] = d/(2^p-1)
    end

    m
end

# Calculate the `w' coefficients in the algorithm by Sweldens and Piessens.
# This is based on Appendix A in the above-mentioned reference.
function wcoefficients(::Type{T}, L, k, M::Int) where {T}
    lambda = T(2)*(k-1)/(T(L)-1)-1

    w = zeros(T,M+1,M+1)

    w[1,1] = 1
    w[2,1] = lambda
    w[2,2] = 1
    if M > 1
        w[3,1] = 2*lambda^2-3
        w[3,2] = 4*lambda
        w[3,3] = 1
    end

    # p+1 is the degree of the next polynomial
    for p in 2:M-1
        w[p+2,1] = w[p+1,2] + 2*lambda*w[p+1,1] - 4*w[p,1]
        w[p+2,2] = 2*w[p+1,1]+w[p+1,3]+2*lambda*w[p+1,2]-4*w[p,2]
        for i in 2:p-1
            w[p+2,i+1] = w[p+1,i]+w[p+1,i+2]+2*lambda*w[p+1,i+1]-4*w[p,i+1]
        end
        w[p+2,p+1] = w[p+1,p] + 2*lambda*w[p+1,p+1]
        w[p+2,p+2] = w[p+1,p+1]
    end

    w
end
