
"Supertype of translates of refinable kernels."
abstract type Refinables{T,EXT} <: BasisTranslates.Translates{T,T,EXT} end

"A basis of translates of a (periodized) refinable function."
struct PeriodicRefinables{T,R<:Refinable{T}} <: Refinables{T,:unitperiodic}
    refinable   ::  R
    n           ::  Int
    # function to be called when evaluating the refinable function
    fun

    function PeriodicRefinables{T,R}(refinable::R, n::Int) where {T,R<:Refinable{T}}
        level = 6
        t, vals = eval_dyadic(coefficients(refinable), level)
        fun = Expansion(RegularBSplines(t, degree=1), vals)
        new(refinable, n, fun)
    end
    function PeriodicRefinables{T,R}(refinable::R, n::Int, fun) where {T,R<:Refinable{T}}
        new(refinable, n, fun)
    end
end

PeriodicRefinables(refinable::R, n::Int) where {T,R<:Refinable{T}} =
    PeriodicRefinables{T,R}(refinable, n)
PeriodicRefinables(refinable::R, n::Int, fun) where {T,R<:Refinable{T}} =
    PeriodicRefinables{T,R}(refinable, n, fun)

Base.size(Φ::PeriodicRefinables) = (Φ.n,)

Base.similar(Φ::PeriodicRefinables{S}, ::Type{S}, n::Int) where {S} =
    PeriodicRefinables(Φ.refinable, n, Φ.fun)

parent_kernel(Φ::PeriodicRefinables) = Φ.refinable

kernel_eval(Φ::PeriodicRefinables, x) = _kernel_eval(Φ.fun, x)
_kernel_eval(fun, x) = fun(x)

function refine(f::Expansion{S,T,B}) where {S,T,B<:PeriodicRefinables}
    n = length(f)
    basis2 = similar(dictionary(f), 2n)
    coef = coefficients(parent_kernel(basis2))
    I = support(coef)
    A = CompactCirculant(datavector(coef), 2n; offset = first(I)+1)
    R = RestrictionArray(2n, 1:2:2n)
    C = R*A
    Expansion(basis2, C'*coefficients(f)*sqrt(S(2)))
end
