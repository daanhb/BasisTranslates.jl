
import BasisFunctions:
    interpolation_grid,
    hasgrid_transform,
    transform_to_grid,
    transform_from_grid

const AllPeriodicTranslates{S,T} =
    Union{PeriodicTranslates{S,T},UnitPeriodicTranslates{S,T}}

period(Φ::UnitPeriodicTranslates{S,T}) where {S,T} = one(S)
support(Φ::UnitPeriodicTranslates{S}) where {S} = UnitInterval{S}()
support(Φ::PeriodicTranslates) = 0..period(Φ)

centers(Φ::UnitPeriodicTranslates{S}) where {S} =
    UnitPeriodicEquispacedGrid{S}(length(Φ))
centers(Φ::PeriodicTranslates) =
    PeriodicEquispacedGrid(support(Φ); length = length(Φ))

hasinterpolationgrid(Φ::AllPeriodicTranslates) = true
interpolation_grid(Φ::AllPeriodicTranslates) = centers(Φ)

"Make a grid with a specified factor of oversampling."
oversampled_grid(Φ::AllPeriodicTranslates; osf::Int) =
    resize(interpolation_grid(Φ), osf*length(Φ))

"""
Without linear scaling the kernel function `φ` is translated to all the centers,
`φ(x-k/n)`. In this case the translates have the same shape regardless of `n`.

With linear scaling the translates are defined by `φ(n*x - k)`. In this case
the support of the translates is inversely proportional to `n`.
"""
linearscaling(Φ::AllPeriodicTranslates) = true

map_to_kernel(Φ::UnitPeriodicTranslates, i) =
    linearscaling(Φ) ? AffineMap{prectype(Φ)}(Φ.n, -(i-1)) : Translation(-center(Φ, i))
function map_to_kernel(Φ::PeriodicTranslates, i)
    c = center(Φ, i)
    if linearscaling(Φ)
        h = step(Φ)
        AffineMap(inv(h), -c/h)
    else
        Translation(-c)
    end
end


function kerneldomain_period(Φ::AllPeriodicTranslates)
    m = map_to_kernel(Φ, 1)
    matrix(m) * period(Φ)
end

function support(Φ::AllPeriodicTranslates, i)
    if hascompactsupport(Φ)
        minv = map_from_kernel(Φ, i)
        a,b = extrema(kernel_support(Φ))
        PeriodicInterval(Interval(minv(a), minv(b)), support(Φ))
    else
        support(Φ)
    end
end

function support_approximate(Φ::AllPeriodicTranslates, i)
    if hascompactsupport_approximate(Φ)
        minv = map_from_kernel(Φ, i)
        a,b = extrema(kernel_support_approximate(Φ))
        PeriodicInterval(Interval(minv(a), minv(b)), support(Φ))
    else
        support(Φ)
    end
end

kernel(Φ::AllPeriodicTranslates) =
    PeriodizedKernel(parent_kernel(Φ), kerneldomain_period(Φ))

kernel_eval(Φ::AllPeriodicTranslates, x) = periodized_kernel_eval(Φ, x)
periodized_kernel_eval(Φ::AllPeriodicTranslates, x) =
    periodized_kernel_eval(parent_kernel(Φ), kerneldomain_period(Φ), x)

kernel_eval_derivative(Φ::AllPeriodicTranslates, x, order) =
    periodized_kernel_eval_derivative(Φ, x, order)
periodized_kernel_eval_derivative(Φ::AllPeriodicTranslates, x, order) =
    periodized_kernel_eval_derivative(parent_kernel(Φ), kerneldomain_period(Φ), x, order)

function undo_periodic_wrapping(Φ::AllPeriodicTranslates, y)
    supp = kernel_support(Φ)
    if y ∈ supp
        y
    else
        a,b = extrema(supp)
        P = kerneldomain_period(Φ)
        y < a ? y+P : y-P
    end
end

function to_kernel_domain(Φ::AllPeriodicTranslates, idx, x)
    m = map_to_kernel(Φ, idx)
    y = undo_periodic_wrapping(Φ, m(x))
    y
end

hasgrid_transform(Φ::UnitPeriodicTranslates, gridspace, grid::UnitPeriodicEquispacedGrid) =
    size(Φ) == size(grid)

transform_from_grid(::Type{T}, src, dest::UnitPeriodicTranslates, grid; options...) where {T} =
    inv(transform_to_grid(T, dest, src, grid; options...))

function transform_to_grid(::Type{T}, src::UnitPeriodicTranslates, dest, grid::AbstractEquispacedGrid; options...) where {T}
    @assert hasgrid_transform(src, dest, grid)
    CirculantOperator{T}(unsafe_eval_element.(Ref(src), 1, grid), src, dest; options...)
end

function evaluation(::Type{T}, Φ::UnitPeriodicTranslates, gb::GridBasis, grid::UnitPeriodicEquispacedGrid;
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
function az_approximate(basis::UnitPeriodicTranslates, f, t1, t2, osf)
    n = length(basis)
    tt = oversampled_grid(basis; osf)

    # Compute the restriction R from the large grid to the subgrid
    I1 = findfirst(tt .>= t1)
    I2 = findlast(tt .<= t2)
    R = RestrictionArray(length(tt), I1:I2)
    tt_int = R*tt

    # Compute the factorization of A using its block-circulant structure
    Acirc = matrix(evaluation(basis, tt))
    F = factorize(Acirc)

    # This leads to an (A,Z) pair on [0,1]
    AZ_A = ArrayOperator(R)*ArrayOperator(F.Π')*ArrayOperator(F.P)*ArrayOperator(F.D)*ArrayOperator(F.F')
    AZ_Zstar = ArrayOperator(F.F)*ArrayOperator(F.Dpinv)*ArrayOperator(F.P')*ArrayOperator(F.Π)*ArrayOperator(R')

    b = f.(tt_int)
    c = az(AZ_A, AZ_Zstar, b)
    F = Expansion(basis, c)
    F
end
