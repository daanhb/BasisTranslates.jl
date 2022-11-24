
import BasisFunctions:
    support,
    isperiodic, period,
    unsafe_eval_element,
    unsafe_eval_element_derivative,
    interpolation_grid,
    hasgrid_transform,
    transform_to_grid,
    transform_from_grid,
    evaluation

"""
Supertype of periodized translates of a single kernel function.

The basis is periodized on the interval [0,1].
"""
abstract type PeriodicTranslates{S,T} <: Translates{S,T} end

isperiodic(Φ::PeriodicTranslates) = true
period(Φ::PeriodicTranslates{S,T}) where {S,T} = one(S)

translates_grid(Φ::PeriodicTranslates{S}) where {S} =
    UnitPeriodicEquispacedGrid{S}(length(Φ))

hasinterpolationgrid(Φ::PeriodicTranslates) = true
interpolation_grid(Φ::PeriodicTranslates) = translates_grid(Φ)

"The map from [0,1] to the domain of the kernel for index `i`."
translate_map(Φ::PeriodicTranslates, i) = AffineMap{prectype(Φ)}(Φ.n, -(i-1))
# By default the map is phi(n*x-k), so that the definition of the
# kernel function may be independent of n.

translate_map_inverse(Φ::PeriodicTranslates, i) = inverse(translate_map(Φ, i))

function kerneldomain_period(Φ::PeriodicTranslates)
    m = translate_map(Φ, 1)
    matrix(m) * period(Φ)
end

function support(Φ::PeriodicTranslates, i)
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

function unsafe_eval_element(Φ::PeriodicTranslates, idx, x)
    m = translate_map(Φ::PeriodicTranslates, idx)
    y = undo_periodic_wrapping(Φ, m(x))
    periodized_kernel_eval(Φ, y)
end

function unsafe_eval_element_derivative(Φ::PeriodicTranslates, idx, x, order)
    m = translate_map(Φ::PeriodicTranslates, idx)
    y = undo_periodic_wrapping(Φ, m(x))
    periodized_kernel_eval_derivative(Φ, order, y)
end


hasgrid_transform(Φ::PeriodicTranslates, gridspace, grid::UnitPeriodicEquispacedGrid) =
    size(Φ) == size(grid)

transform_from_grid(::Type{T}, src, dest::PeriodicTranslates, grid; options...) where {T} =
    inv(transform_to_grid(T, dest, src, grid; options...))

function transform_to_grid(::Type{T}, src::PeriodicTranslates, dest, grid::AbstractEquispacedGrid; options...) where {T}
    @assert hasgrid_transform(src, dest, grid)
    CirculantOperator{T}(unsafe_eval_element.(Ref(src), 1, grid), src, dest; options...)
end

function evaluation(::Type{T}, Φ::PeriodicTranslates, gb::GridBasis, grid::UnitPeriodicEquispacedGrid;
            options...) where {T}
    n_grid = length(grid)
    n_Φ = length(Φ)
    s, rem = divrem(n_grid, n_Φ)
    if rem == 0
        h_small = one(real(T))/n_grid
        h_large = one(real(T))/n_Φ
        circulant_arrays = [Circulant{T}(unsafe_eval_element.(Ref(Φ), 1, range(k*h_small, length=n_Φ, step=h_large))) for k in 0:s-1]
        M = MultiRowCirculant(circulant_arrays)
        ArrayOperator{T}(M, Φ, gb)
    else
        BasisFunctions.dense_evaluation(T, Φ, gb; options...)
    end
end
