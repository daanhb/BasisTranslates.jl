"Supertype of radial basis functions."
abstract type RBF{T}
end

"A smooth radial basis function, typically used with a shape parameter `ε`."
abstract type SmoothRBF{T} <: RBF{T}
end

hascompactsupport(rbf::SmoothRBF) = false

"A piecewise smooth radial basis function."
abstract type PiecewiseSmoothRBF{T} <: RBF{T}
end

"The Multiquadric RBF function."
struct Multiquadric{T} <: SmoothRBF{T}
end
const MQ{T} = Multiquadric{T}
MQ() = MQ{Float64}()

"The Inverse Multiquadric RBF function."
struct InverseMultiquadric{T} <: SmoothRBF{T}
end
const IMQ{T} = InverseMultiquadric{T}
IMQ() = IMQ{Float64}()

"The Inverse Quadratic RBF function."
struct InverseQuadratic{T} <: SmoothRBF{T}
end
const IQ{T} = InverseQuadratic{T}
IQ() = IQ{Float64}()

struct Gaussian{T} <: SmoothRBF{T}
end
const GA{T} = Gaussian{T}
GA() = GA{Float64}()

approximate_support(rbf::Gaussian{T}, threshold = eps(T)) where {T} =
    sqrt(-log(threshold))


"The polyharmonic spline RBF function."
struct PolyharmonicSpline{T} <: PiecewiseSmoothRBF{T}
    p   ::  Int
end
const PHS{T} = PolyharmonicSpline{T}

PolyharmonicSpline(p::Int) = PolyharmonicSpline{Float64}(p)

(rbf::PolyharmonicSpline)(r) = rbf_phs(r, r.p)

similar(r::PolyharmonicSpline{T}, p) where {T} = PolyharmonicSpline{T}(p)


"The thin plate spline RBF function."
struct ThinPlateSpline{T} <: PiecewiseSmoothRBF{T}
    p   ::  Int
end
const TPS{T} = ThinPlateSpline{T}

ThinPlateSpline(p::Int) = ThinPlateSpline{Float64}(p)

(rbf::ThinPlateSpline)(r) = rbf_tps(r, rbf.p)

similar(r::ThinPlateSpline{T}, p) where {T} = ThinPlateSpline{T}(p)
