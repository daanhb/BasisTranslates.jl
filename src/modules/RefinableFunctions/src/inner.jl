
"A matrix with the convolution of two sequences on its rows, shifted by two on consecutive rows."
function transfer_matrix(s1::CompactSequence, s2::CompactSequence)
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

# "Return the integrals `φ1(x)*φ2(x-k)` for all k."
# function unscaled_productintegral_moments(s1::CompactSequence, s2::CompactSequence)
#     A,I = transfer_matrix(s1, s2)
#     vals = eigenvector_with_eigenvalue_1(A)
#     # TODO: figure out where the range of the moments is
#     VectorSequence(vals, first(I)+2-datalength(s1):last(I)-datalength(s1))
# end
