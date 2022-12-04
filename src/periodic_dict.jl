
import BasisFunctions:
    interpolation_grid,
    hasgrid_transform,
    transform_to_grid,
    transform_from_grid

"""
Supertype of periodized translates of a single kernel function.
"""
abstract type PeriodicTranslates{S,T} <: Translates{S,T} end

isperiodic(Φ::PeriodicTranslates) = true
period(Φ::PeriodicTranslates{S,T}) where {S,T} = one(S)

support(Φ::PeriodicTranslates{S}) where {S} = UnitInterval{S}()

centers(Φ::PeriodicTranslates{S}) where {S} =
    UnitPeriodicEquispacedGrid{S}(length(Φ))

hasinterpolationgrid(Φ::PeriodicTranslates) = true
interpolation_grid(Φ::PeriodicTranslates) = centers(Φ)

oversampling_grid(Φ::PeriodicTranslates, s::Int) =
    resize(interpolation_grid(Φ), s*length(Φ))

"""
Without linear scaling the kernel function `φ` is translated to all the centers,
`φ(x-k/n)`. In this case the translates have the same shape regardless of `n`.

With linear scaling the translates are defined by `φ(n*x - k)`. In this case
the support of the translates is inversely proportional to `n`.
"""
linearscaling(Φ::PeriodicTranslates) = true

map_to_kernel(Φ::PeriodicTranslates, i) =
    linearscaling(Φ) ? AffineMap{prectype(Φ)}(Φ.n, -(i-1)) : Translation(-(i-one(prectype(Φ)))/length(Φ))

function kerneldomain_period(Φ::PeriodicTranslates)
    m = map_to_kernel(Φ, 1)
    matrix(m) * period(Φ)
end

function support(Φ::PeriodicTranslates, i)
    if hascompactsupport(Φ)
        minv = map_from_kernel(Φ, i)
        a,b = extrema(kernel_support(Φ))
        PeriodicInterval(Interval(minv(a), minv(b)), support(Φ))
    else
        support(Φ)
    end
end

function support_approximate(Φ::PeriodicTranslates, i)
    if hascompactsupport_approximate(Φ)
        minv = map_from_kernel(Φ, i)
        a,b = extrema(kernel_support_approximate(Φ))
        PeriodicInterval(Interval(minv(a), minv(b)), support(Φ))
    else
        support(Φ)
    end
end

kernel(Φ::PeriodicTranslates) =
    PeriodizedKernel(parent_kernel(Φ), kerneldomain_period(Φ))

kernel_eval(Φ::PeriodicTranslates, x) = periodized_kernel_eval(Φ, x)
periodized_kernel_eval(Φ::PeriodicTranslates, x) =
    periodized_kernel_eval(parent_kernel(Φ), kerneldomain_period(Φ), x)

kernel_eval_derivative(Φ::PeriodicTranslates, x, order) =
    periodized_kernel_eval_derivative(Φ, x, order)
periodized_kernel_eval_derivative(Φ::PeriodicTranslates, x, order) =
    periodized_kernel_eval_derivative(parent_kernel(Φ), kerneldomain_period(Φ), x, order)

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

function to_kernel_domain(Φ::PeriodicTranslates, idx, x)
    m = map_to_kernel(Φ, idx)
    y = undo_periodic_wrapping(Φ, m(x))
    y
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

"""
Approximate the given function `f` on the interval `[t1,t2] ⊂ [0,1]` using `n`
centers on `[0,1]` and with a factor `osf` of oversampling.
"""
function az_approximate(basis::PeriodicTranslates, f, t1, t2, osf)
    n = length(basis)
    tt = oversampling_grid(basis, osf)

    # Compute the restriction R from the large grid to the subgrid
    I1 = findfirst(tt .>= t1)
    I2 = findlast(tt .<= t2)
    R = ContiguousRows(length(tt), I1:I2)
    tt_int = R*tt

    # Compute the factorization of A using its block-circulant structure
    Acirc = matrix(evaluation(basis, tt))
    F = factorize(Acirc)

    # This leads to an (A,Z) pair on [0,1]
    AZ_A = R*(F.Π'*(F.P*(F.D*LinearMap(F.F'))))
    AZ_Zstar = F.F*(F.Dpinv*(F.P'*(F.Π*LinearMap(R'))))

    b = f.(tt_int)
    c = az(AZ_A, AZ_Zstar, b)
    F = Expansion(basis, c)

    F
end
