
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


"""
Supertype of periodized translates of a single kernel function.

The basis is periodized on the interval [0,1].
"""
abstract type PeriodicTranslates{S,T} <: Translates{S,T} end

BasisFunctions.isperiodic(Φ::PeriodicTranslates) = true
BasisFunctions.period(Φ::PeriodicTranslates{S,T}) where {S,T} = one(S)

translates_grid(Φ::PeriodicTranslates{S}) where {S} =
    UnitPeriodicEquispacedGrid{S}(length(Φ))

"The map from [0,1] to the domain of the kernel for index `i`."
translate_map(Φ::PeriodicTranslates, i) = AffineMap{prectype(Φ)}(Φ.n, -(i-1))
# By default the map is phi(n*x-k), so that the definition of the
# kernel function may be independent of n.

translate_map_inverse(Φ::PeriodicTranslates, i) = inverse(translate_map(Φ, i))

function kerneldomain_period(Φ::PeriodicTranslates)
    m = translate_map(Φ, 1)
    matrix(m) * period(Φ)
end

function BasisFunctions.support(Φ::PeriodicTranslates, i)
    if hascompactsupport(Φ)
        minv = translate_map_inverse(Φ, i)
        a,b = extrema(kernel_support(Φ))
        PeriodicInterval(Interval(minv(a), minv(b)), support(Φ))
    else
        support(Φ)
    end
end

function support_approximate(Φ::PeriodicTranslates, i)
    if hascompactsupport_approximate(Φ)
        minv = translate_map_inverse(Φ, i)
        a,b = extrema(kernel_support_approximate(Φ))
        PeriodicInterval(Interval(minv(a), minv(b)), support(Φ))
    else
        support(Φ)
    end
end

"""
    k_left, k_right = nb_overlapping_kernels(Φ::PeriodicTranslates)

The number of shifted kernels that overlap, due to periodicity. The functions
`phi(x + k*P)` may overlap with `phi(x)` for all `k` in `1:k_left`, and similarly
for `phi(x - k*P)` and `1:k_right`.
"""
function nb_overlapping_kernels(Φ::PeriodicTranslates)
    a,b = extrema(kernel_support_approximate(Φ))
    P = kerneldomain_period(Φ)
    k = floor(Int, (b-a)/P)
    k, k
end

function nb_overlapping_kernels(Φ::PeriodicTranslates, x)
    a,b = extrema(kernel_support_approximate(Φ))
    P = kerneldomain_period(Φ)
    k_left = floor(Int, (b-x)/P)
    k_right = floor(Int, (x-a)/P)
    k_left, k_right
end

"""
Evaluate the periodized kernel of the basis of translates.

It is defined as the infinite sum `sum_k( phi(x + k*P))`, where `P`
is the period in the kernel domain.
"""
function periodized_kernel_eval(Φ::PeriodicTranslates, x)
    k_left, k_right = nb_overlapping_kernels(Φ, x)
    P = kerneldomain_period(Φ)
    z = kernel_eval(Φ, x)
    for k in 1:k_right
        z1 = kernel_eval(Φ, x-k*P)
        z += z1
    end
    for k in 1:k_left
        z1 = kernel_eval(Φ, x+k*P)
        z += z1
    end
    z
end

"Evaluate the derivative of the periodized kernel."
function periodized_kernel_eval_derivative(Φ::PeriodicTranslates, order, x)
    k_left, k_right = nb_overlapping_kernels(Φ, x)
    P = kerneldomain_period(Φ)
    z = kernel_eval_derivative(Φ, order, x)
    for k in 1:k_right
        z1 = kernel_eval_derivative(Φ, order, x-k*P)
        z += z1
    end
    for k in 1:k_left
        z1 = kernel_eval_derivative(Φ, order, x+k*P)
        z += z1
    end
    z
end

function undo_periodic_wrapping(Φ::PeriodicTranslates, y)
    supp = kernel_support(Φ)
    if y ∈ supp
        y
    else
        a,b = extrema(supp)
        P = kerneldomain_period(Φ)
        y < a ? y+P : y-P
    end
end

function BasisFunctions.unsafe_eval_element(Φ::PeriodicTranslates, idx, x)
    m = translate_map(Φ::PeriodicTranslates, idx)
    y = undo_periodic_wrapping(Φ, m(x))
    periodized_kernel_eval(Φ, y)
end

function BasisFunctions.unsafe_eval_element_derivative(Φ::PeriodicTranslates, idx, x, order)
    m = translate_map(Φ::PeriodicTranslates, idx)
    y = undo_periodic_wrapping(Φ, m(x))
    periodized_kernel_eval_derivative(Φ, order, y)
end
