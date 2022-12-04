
import BasisFunctions:
    support,
    unsafe_eval_element,
    unsafe_eval_element_derivative,
    evaluation

"Supertype of a basis of translates."
abstract type Translates{S,T} <: Dictionary{S,T} end

"Return the centers of the basis of translates."
function centers() end

"Return the center of the basis function with index `i`."
center(Φ::Translates, i) = centers(Φ)[i]

Base.size(Φ::Translates) = size(centers(Φ))

# The natural index is that of the grid of centers
BasisFunctions.ordering(Φ::Translates) = eachindex(centers(Φ))

"Map from the domain of the basis function to the kernel domain."
map_to_kernel(Φ::Translates, i) = Translation(-center(Φ, i))
"Inverse of `map_to_kernel`."
map_from_kernel(Φ::Translates, i) = inverse(map_to_kernel(Φ, i))

"Map the given point `x` from the basis domain to the kernel domain."
to_kernel_domain(Φ::Translates, idx, x) = map_to_kernel(Φ, idx)(x)
"Map the given point `y` from the kernel domain to the basis domain."
from_kernel_domain(Φ::Translates, idx, y) = map_from_kernel(Φ, idx)(y)

unsafe_eval_element(Φ::Translates, idx, x) =
    kernel_eval(Φ, to_kernel_domain(Φ, idx, x))

# In the computation of the derivative, we have to take into account the
# jacobian of the map from the basis domain to the kernel domain
unsafe_eval_element_derivative(Φ::Translates, idx, x, order) =
    _eval_deriv(map_to_kernel(Φ, idx), Φ, idx, x, order)
_eval_deriv(m::Translation, Φ, idx, x, order) =
    kernel_eval_derivative(Φ, to_kernel_domain(Φ, idx, x), order)
function _eval_deriv(m, Φ, idx, x, order)
    @assert isaffine(m)
    jacobian(m, x)^order * kernel_eval_derivative(Φ, to_kernel_domain(Φ, idx, x), order)
end

function support(Φ::Translates{S}) where {S}
    if hascompactsupport(Φ)
        a,b = extrema(kernel_support(Φ))
        minv1 = map_from_kernel(Φ, 1)
        minv2 = map_from_kernel(Φ, length(Φ))
        minv1(a)..minv2(b)
    else
        FullSpace{S}()
    end
end

function support(Φ::Translates, i)
    if hascompactsupport(Φ)
        minv = map_from_kernel(Φ, i)
        a,b = extrema(kernel_support(Φ))
        (minv(a)..minv(b)) ∩ support(Φ)
    else
        support(Φ)
    end
end

BasisFunctions.hasmeasure(Φ::Translates) = true
BasisFunctions.measure(Φ::Translates) = lebesguemeasure(support(Φ))

"Evaluate the kernel of the basis of translates."
kernel_eval(Φ::Translates{S}, x) where {S} = kernel_eval(kernel(Φ), S(x))

"Evaluate the derivative of the kernel of the basis of translates."
kernel_eval_derivative(Φ::Translates{S}, x, order) where {S} =
    kernel_eval_derivative(kernel(Φ), S(x), order)

"The support of the kernel of the basis of translates."
kernel_support(Φ::Translates) = kernel_support(kernel(Φ))

"The approximate support of the kernel based on thresholding."
kernel_support_approximate(Φ::Translates, threshold = eps(prectype(Φ))) =
    kernel_support_approximate(kernel(Φ), threshold)

"Do the basis functions have compact support?"
hascompactsupport(Φ::Translates) = iscompact(kernel_support(Φ))

"Do the basis functions approximately have compact support?"
hascompactsupport_approximate(Φ::Translates, threshold...) =
    iscompact(kernel_support_approximate(Φ, threshold...))

kerneldomain_center(Φ::Translates, i) =
    to_kernel_domain(Φ, center(Φ, i))
