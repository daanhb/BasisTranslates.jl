
lightning_uniform_poles(σ::Number, n) = lightning_uniform_poles(float(σ), n)
lightning_uniform_poles(σ::T, n) where {T <: AbstractFloat} =
    -exp.(-σ*(0:n-1)./sqrt(T(n)))

lightning_tapered_poles(σ::Number, n) = lightning_tapered_poles(float(σ), n)
lightning_tapered_poles(σ::T, n) where {T <: AbstractFloat} =
    -exp.(σ*(sqrt.(T(n):-1:1) .- sqrt(T(n))))

function lightning_samples(poles::AbstractVector{T}, osf::Int) where {T}
    n = length(poles)
    m = osf*n
    samples = zeros(T, m)
    pts = [0; sort(-poles)]
    for i in 1:n
        samples[(i-1)*osf+1:i*osf] = range(pts[i], stop=pts[i+1], length=osf+2)[2:end-1]
    end
    samples
end

lightning_normalization(Φ::BasisTranslates.Translates) =
    DiagonalOperator(Φ, -centers(Φ))


function poly_lightning_approx(f, n, ncheb_uni, ncheb_tap, osf, kernel, sigma_uni = sqrt(2)*pi, sigma_tap = 2*sqrt(2)*pi; C = 1)
    poles_uni = lightning_uniform_poles(sigma_uni, n)
    poles_tap = lightning_tapered_poles(sigma_tap, n)

    basis_uni = KernelTranslates(kernel, C*poles_uni)
    basis_tap = KernelTranslates(kernel, C*poles_tap)
    if kernel isa SingKernel
        nbasis_uni = DiagonalOperator(basis_uni, [1/bf(0) for bf in basis_uni])*basis_uni
        nbasis_tap = DiagonalOperator(basis_tap, [1/bf(0) for bf in basis_tap])*basis_tap
    else
        nbasis_uni = basis_uni
        nbasis_tap = basis_tap
    end

    pts_uni = lightning_samples(poles_uni, osf)
    pts_tap = lightning_samples(poles_tap, osf)
    if ncheb_uni > 0
        chebpts_uni = (ChebyshevNodes(ncheb_uni*osf) .+ 1)/2
        chebpts_tap = (ChebyshevNodes(ncheb_tap*osf) .+ 1)/2
        pts_uni = [pts_uni; chebpts_uni]
        pts_tap = [pts_tap; chebpts_tap]
        cheb_uni = ChebyshevT(ncheb_uni) → UnitInterval()
        cheb_tap = ChebyshevT(ncheb_tap) → UnitInterval()
        lightning_basis_uni = nbasis_uni ⊕ cheb_uni
        lightning_basis_tap = nbasis_tap ⊕ cheb_tap
    else
        lightning_basis_uni = nbasis_uni
        lightning_basis_tap = nbasis_tap
    end
    A_uni = evaluation_matrix(lightning_basis_uni, pts_uni)
    b_uni = f.(pts_uni)
    c_uni = A_uni\b_uni
    c2 = zeros(lightning_basis_uni)
    c2 .= c_uni
    F_uni = Expansion(lightning_basis_uni, c2)
    A_tap = evaluation_matrix(lightning_basis_tap, pts_tap)
    b_tap = f.(pts_tap)
    c_tap = A_tap\b_tap
    c2 = zeros(lightning_basis_tap)
    c2 .= c_tap
    F_tap = Expansion(lightning_basis_tap, c2)
    F_uni, F_tap, poles_uni, poles_tap
end

function maxerror(f, F, osf, poles, ncheb)
    n = length(poles)
    pts_poles = lightning_samples(poles, osf)
    if ncheb > 0
        pts_cheb = (ChebyshevNodes(ncheb*osf) .+ 1)/2
        pts = [pts_poles; pts_cheb]
    else
        pts = pts_poles
    end
    maximum(abs.(f.(pts)-F.(pts)))
end

struct SingKernel
    α :: Float64
end
kernel_eval(k::SingKernel, x) = x^(k.α)
