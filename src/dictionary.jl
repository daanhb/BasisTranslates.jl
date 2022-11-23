
"Supertype of a basis of translates."
abstract type Translates{S,T} <: Dictionary{S,T} end

"Return the grid of the basis of translates."
function translates_grid() end

translate_center(Φ::Translates, i) = translates_grid(Φ)[i]

# The natural index is that of the grid of centers
BasisFunctions.ordering(Φ::Translates) = eachindex(translates_grid(Φ))

function BasisFunctions.unsafe_eval_element(Φ::Translates{S,T}, idx, x) where {S,T}
    m = translate_map(Φ, idx)
    kernel_eval(Φ, m(x))
end

BasisFunctions.support(Φ::Translates) = coverdomain(translates_grid(Φ))
BasisFunctions.hasmeasure(Φ::Translates) = true
BasisFunctions.measure(Φ::Translates) = lebesguemeasure(support(Φ))

"Evaluate the kernel of the basis of translates."
function kernel_eval() end
"Evaluate the derivative of the kernel of the basis of translates."
function kernel_eval_derivative() end
"The support of the kernel of the basis of translates."
kernel_support(Φ::Translates) = FullSpace{domaintype(Φ)}()

"The approximate support of the kernel based on thresholding."
kernel_support_approximate(Φ::Translates, threshold = eps(prectype(Φ))) =
    kernel_support(Φ)

"Do the basis functions have compact support?"
hascompactsupport(Φ::Translates) = kernel_support(Φ) isa Interval

"Do the basis functions approximately have compact support?"
hascompactsupport_approximate(Φ::Translates, threshold...) =
    kernel_support_approximate(Φ, threshold...) isa Interval

kerneldomain_center(Φ::Translates, i = 1) = translate_map(Φ, i)(0)
