
"""
Construct the array representing periodic convolution of one compact
sequence with another one.
"""
function periodic_filter_array(s::CompactSequence, n::Int)
    srev = reverse(s)
    vec = datavector(srev)
    I = support(srev)
    CompactCirculant(vec, n, offset = first(I)+1)
end

function periodic_downsampled_filter_array(s::CompactSequence, n::Int)
    R = RestrictionArray(n, 1:2:n)
    A = periodic_filter_array(s, n)
    R*A
end

periodic_upsampling_filter_array(s::CompactSequence, n::Int) =
    adjoint(periodic_downsampled_filter_array(reverse(s), n))
