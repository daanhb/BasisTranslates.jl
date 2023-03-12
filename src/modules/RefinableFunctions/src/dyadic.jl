
"Return the eigenvector of the matrix `A` with eigenvalue `1`."
function eigenvector_with_eigenvalue_1(A)
    e, v = eigen(A)
    # There should be an eigenvalue close to 1
    @assert any(abs.(e .- 1) .< 1e-10)
    I = findfirst(abs.(e .- 1) .< 1e-10)
    v[:,I]
end

"""
Use an eigenvalue problem to find the exact values of the refinable function
at integer points.

The result is normalized by ensuring the sum of values equals the given moment.
"""
function values_at_integers(s::CompactSequence, moment = 1)
    I = support(s)
    i1 = first(I)
    L = length(I)
    A = zeros(eltype(s),L,L)
    for l in 0:L-1
        for k in 0:L-1
            A[k+1,l+1] = sqrt(2)*s[i1-l+2k]
        end
    end
    vals = eigenvector_with_eigenvalue_1(A)
    vals/sum(vals)*moment
end

"""
Use the two-scale relation to find the function values at the next level of
dyadic numbers.
"""
function refine_dyadic_values(s::CompactSequence, vals::Vector)
    I = support(s)
    # we shift the sequence below by i1 to make it causal
    i1 = first(I)
    L = length(I)
    level = round(Int, (length(vals)-1)/(L-1))-1
    newvals = similar(vals, 2*length(vals)-1)
    newvals[1:2:end] = vals
    # the values in vals satisfy vals[j] = φ(2^{-level}*(j-1))
    # the values in newvals satisfy newvals[j] = φ(2^{-level-1}*(j-1))
    for k in 2:2:length(newvals)-1
        z = zero(eltype(vals))
        for i in 0:L-1
            val_idx = k - (level+1)*i
            if 1 <= val_idx <= length(vals)
                z = z + sqrt(2)*s[i1+i]*vals[val_idx]
            end
        end
        newvals[k] = z
    end
    newvals
end

"Evaluate a refinable function at dyadic numbers up to the given level."
function eval_dyadic(s::CompactSequence, level::Int, moment = 1)
    vals = values_at_integers(s, moment)
    for k in 1:level
        vals = refine_dyadic_values(s, vals)
    end
    I = support(s)
    T = eltype(s)
    t = first(I):T(2)^(-level):last(I)
    t, vals
end
