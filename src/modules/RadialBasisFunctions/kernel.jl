"Supertype of radial basis functions."
abstract type RBF
end

function kernel_eval_derivative(rbf::RBF, order::Int, r)
    if order == 0
        kernel_eval(rbf, r)
    elseif order == 1
        kernel_eval_diff1(rbf, r)
    elseif order == 2
        kernel_eval_diff2(rbf, r)
    elseif order == 3
        kernel_eval_diff3(rbf, r)
    elseif order == 4
        kernel_eval_diff4(rbf, r)
    else
        error("High order differentiation not implemented.")
    end
end

"A smooth radial basis function, typically used with a shape parameter `ε`."
abstract type SmoothRBF <: RBF
end

hascompactsupport(rbf::SmoothRBF) = false

"A piecewise smooth radial basis function."
abstract type PiecewiseSmoothRBF <: RBF
end

"The Multiquadric RBF function."
struct Multiquadric <: SmoothRBF
end
const MQ = Multiquadric

"The Inverse Multiquadric RBF function."
struct InverseMultiquadric <: SmoothRBF
end
const IMQ = InverseMultiquadric

"The Inverse Quadratic RBF function."
struct InverseQuadratic <: SmoothRBF
end
const IQ = InverseQuadratic

kernel_eval(rbf::IQ, r) = rbf_iq(r)

struct Gaussian <: SmoothRBF
end
const GA = Gaussian

function support_approximate(rbf::Gaussian, threshold = eps(Float64))
    S = sqrt(-log(threshold))
    -S..S
end

kernel_eval(rbf::GA, r) = rbf_ga(r)
kernel_eval_diff1(rbf::GA, r) = rbf_ga_diff1(r)
kernel_eval_diff2(rbf::GA, r) = rbf_ga_diff2(r)


"The polyharmonic spline RBF function."
struct PolyharmonicSpline <: PiecewiseSmoothRBF
    p   ::  Int
end
const PHS = PolyharmonicSpline

kernel_eval(rbf::PolyharmonicSpline, r) = rbf_phs(r, r.p)


"The thin plate spline RBF function."
struct ThinPlateSpline <: PiecewiseSmoothRBF
    p   ::  Int
end
const TPS = ThinPlateSpline

kernel_eval(rbf::ThinPlateSpline, r) = rbf_tps(r, rbf.p)
