
"""
A dictionary for the (α,β) formulation of a rational function in
barycentric form:

`sum(α[k]/(z-z[k])) / sum(β[k]/(z-z[k]))`

The elements of the dictionary correspond to the individual terms
`1/(z-z[j])` in the numerator, divided by the full denominator. The
coefficients of an expansion therefore determine the `α` vector.

The `β` values are commonly called the weights.
"""
struct BarycentricRational{T} <: BasisFunctions.Dictionary{T,T}
    z   ::  Vector{T}   # support points
    w   ::  Vector{T}   # weights

    function BarycentricRational{T}(z, w) where T
        @assert length(z) == length(w)
        new(z, w)
    end
end

function BarycentricRational(z, w)
    T = promote_type(eltype(z),eltype(w))
    BarycentricRational{T}(z, w)
end

BasisFunctions.size(Φ::BarycentricRational) = size(Φ.w)

BasisFunctions.support(Φ::BarycentricRational{T}) where T =
    DomainSets.FullSpace{T}()

function BasisFunctions.similar(Φ::BarycentricRational, ::Type{T}, n::Int) where T
    @assert n == length(Φ)
    BarycentricRational{T}(Φ.z, Φ.w)
end
Base.complex(Φ::BarycentricRational{T}) where T = similar(Φ, complex(T))

BasisFunctions.interpolation_grid(Φ::BarycentricRational) = Φ.z

"""
Construct the interpolation operator of the barycentric form.

This is a diagonal operator, since the weights of the numerator are simply those
of the denominator times the function values.
"""
function interpolation_operator(Φ::BarycentricRational)
    pts = Φ.z
    DiagonalOperator(GridBasis(pts), Φ, Φ.w)
end

function BasisFunctions.interpolation(::Type{T}, Φ::BarycentricRational, gb::GridBasis) where T
    if grid(gb) == Φ.z
        interpolation_operator(similar(Φ, T))
    else
        BasisFunctions.default_interpolation(T, Φ, gb)
    end
end

"""
Evaluate a rational function in `(α,β)` form with the given support points
'z' in the point `x`.
"""
function barycentric_eval(z, α, β, x)
    if isinf(x)
        r = sum(α) / sum(β)
    else
        r = pf_eval(z, α, x) / pf_eval(z, β, x)
        if isnan(r)
            I = findfirst(x .== z)
            r = α[I]/β[I]
        end
    end
    r
end

"Evaluate one term of a rational function in `(α,β)` form."
function barycentric_eval_term(z, w, i, x)
    num = 1 / (x-z[i])
    den = pf_eval(z, w, x)
    num/den
end

BasisFunctions.unsafe_eval_element(Φ::BarycentricRational, i, x) =
    barycentric_eval_term(Φ.z, Φ.w, i, x)

BasisFunctions.unsafe_eval_expansion(Φ::BarycentricRational, coefficients, x) =
    barycentric_eval(Φ.z, coefficients, Φ.w, x)


"Find the roots of a rational function in partial fractions form."
function pf_roots(z, w)
    m = length(w)
    T = promote_type(eltype(z),eltype(w))
    B = zeros(T,m+1,m+1)
    for i in 2:m+1
        B[i,i] = 1
    end
    E = [0 transpose(w); ones(T, m, 1) Diagonal(z)]
    if T == Float64
        roots,~ = eigen(E, B)
        roots = filter(!isinf, roots)
    else
        EB = inv(E)*B
        pp = eigvals(EB)
        large_enough = abs.(pp) .> 1e-10
        roots = 1 ./ pp[large_enough]
    end
    roots
end

"""
Compute the poles, residues and zeros of a rational function in
barycentric `(α, β)` form.
"""
function barycentric_poles_residues_zeros(z, α, β)
    poles = pf_roots(z, β)  # roots of the denominator
    zeros = pf_roots(z, α)  # roots of the numerator

    # we obtain the residue by differentiation
    den_diff = t -> sum(-β[k]/(t-z[k])^2 for k in 1:length(z))
    residues = pf_eval.(Ref(z), Ref(α), poles) ./ den_diff.(poles)

    poles, residues, zeros
end

"""
Compute the poles of a rational function in barycentric form with the
given weights.
"""
barycentric_poles(z, w) = pf_roots(z, w)

poles(Φ::BarycentricRational) = barycentric_poles(Φ.z, Φ.w)
