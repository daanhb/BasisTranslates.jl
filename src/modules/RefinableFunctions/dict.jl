
"A basis of translates of a (periodized) refinable function."
struct PeriodicRefinables{T,R} <: BasisTranslates.PeriodicTranslates{T,T}
    refinable   ::  R
    n           ::  Int
end
