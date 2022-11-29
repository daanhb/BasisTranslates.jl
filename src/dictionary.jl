
import BasisFunctions:
    support,
    unsafe_eval_element,
    unsafe_eval_element_derivative,
    evaluation

"Supertype of a basis of translates."
abstract type Translates{S,T} <: Dictionary{S,T} end

"Return the grid of the basis of translates."
function translates_grid() end

translate_center(Φ::Translates, i) = translates_grid(Φ)[i]

Base.size(Φ::Translates) = (length(translates_grid(Φ)),)

# The natural index is that of the grid of centers
BasisFunctions.ordering(Φ::Translates) = eachindex(translates_grid(Φ))

"Map from the domain of the basis to the kernel domain for index `i`."
translate_map(Φ::Translates, i) = Translation(-translate_center(Φ, i))

to_kernel_domain(Φ::Translates, idx, x) = x - translate_center(Φ, idx)

unsafe_eval_element(Φ::Translates, idx, x) =
    kernel_eval(Φ, to_kernel_domain(Φ, idx, x))

function unsafe_eval_element_derivative(Φ::Translates, idx, x, order)
    # default implementation is only valid for translations
    m = translate_map(Φ, idx)::Translation
    kernel_eval_derivative(Φ, m(x), order)
end

support(Φ::Translates) = coverdomain(translates_grid(Φ))

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

kerneldomain_translate_center(Φ::Translates, i) =
    to_kernel_domain(Φ, translate_center(Φ, i))
