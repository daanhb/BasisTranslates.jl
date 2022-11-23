
using BasisTranslates.BSplines

struct BSpline{T}
    order   ::  Int
end
(f::BSpline{T})(x) where {T} = eval_bspline(f.order, x, T)

struct PeriodicBSpline{T}
    order   ::  Int
    period  ::  T
end
PeriodicBSpline{T}(order) where {T} = PeriodicBSpline{T}(order, 1)
(f::PeriodicBSpline{T})(x) where {T} =
    eval_periodic_bspline(f.order, x, f.period, T)

struct CenteredBSpline{T}
    order   ::  Int
end
(f::CenteredBSpline{T})(x) where {T} =
    eval_centered_bspline(f.order, x, T)

struct PeriodicCenteredBSpline{T}
    order   ::  Int
    period  ::  T
end
PeriodicCenteredBSpline{T}(order) where {T} =
    PeriodicCenteredBSpline{T}(order, 1)

(f::PeriodicCenteredBSpline{T})(x) where {T} =
    eval_periodic_centered_bspline(f.order, x, f.period, T)


function elementary_props_of_splines_test(T)
    T = real(T)
    tol = sqrt(eps(T))
    S = 20
    for N in 1:10
        f = BSpline{T}(N-1)
        # Integral should be 1
        if !(T <: BigFloat)
            I,e = QuadGK.quadgk(f, 0, N, rtol = tol)
            @test I≈T(1)
        end
        # Infinite summation of shifted versions is 1
        xx = LinRange(T(N-1), T(N), S)[1:end-1]
        # xx = range(T(N-1), step=T(1)/T(S), length=S)
        g = zeros(T,length(xx))
        for k in 0:N-1
            g += map(x->f(x-k), xx)
        end
        @test g ≈ ones(T,length(g))  # (norm(g-1) < tol)
        # Two scale relation
        x = LinRange(T(-1), T(N+1), S)
        # x = range(T(-1), step=T(1)/T(S), length=S)
        g = zeros(T,length(x))
        for k in 0:N
            g += T(binomial(N,k))*map(x->f(2x-k), x)
        end
        g *= T(2)^(-N+1)
        G = map(f, x)
        @test g ≈ G
    end
end

function periodic_bsplines_test(T)
    for N in 0:4
        period = T(N+1)
        x = LinRange(T(0),period,10)[1:end-1]
        PeriodicBSpline{T}(N, period).(x) ≈ BSpline{T}(N).(x)
        for k in -2:2
            @test PeriodicBSpline{T}(N, period).(x) ≈ PeriodicBSpline{T}(N,period).(x .+ k*period)
        end
        period = T(N+1)/T(3)
        x = LinRange(T(0),period,10)[1:end-1]
        for k in -2:2
            @test PeriodicBSpline{T}(N, period).(x) ≈ PeriodicBSpline{T}(N,period).(x .+ k*period)
        end
    end
end

function centered_bsplines_test(T)
    K = 20
    for N in 0:4
        x = LinRange(T(0)+eps(T), T(1), K)
        @test PeriodicCenteredBSpline{T}(N, 10(N+1)).(x) ≈ PeriodicCenteredBSpline{T}(N, 10(N+1)).(-x)
        @test CenteredBSpline{T}(N).(x) ≈ CenteredBSpline{T}(N).(-x)
    end
end

using BasisTranslates.BSplines:
    squared_spline_integral, shifted_spline_integral

function test_spline_integrals()
    T = Float64
    for N in 0:8
        @test squared_spline_integral(N) ≈ shifted_spline_integral(N,0)
        @test squared_spline_integral(N, Float64) ≈ shifted_spline_integral(N,0, Float64)
        @test squared_spline_integral(N) ≈ shifted_spline_integral(N,0, Float64)
        for t in 0:4
            f = BSpline{T}(N)
            i = QuadGK.quadgk(x->f(x)*f(x+t),(-N-2:N+2)..., rtol=sqrt(eps(T)))[1]
            @test abs(i-shifted_spline_integral(N,t))<10*sqrt(eps(T))
            @test abs(i-shifted_spline_integral(N,t,Float64))<10*sqrt(eps(T))
            @test abs(i-shifted_spline_integral(N,t,BigFloat))<10*sqrt(eps(T))
        end
    end
end


function derivative_test(T)
    t = LinRange(-10,10,100)
    t2 = (t[1:end-1]+t[2:end])/2
    for i in 2:10
        y1 = eval_bspline.(i, t, T)
        y2 = eval_bspline_derivative.(i, 1, t2, T)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end
    for i in 2:10
        y1 = eval_periodic_bspline.(i, t, 5, T)
        y2 = eval_periodic_bspline_derivative.(i, 1, t2, 5, T)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end

    for t in [LinRange(0,1,100), LinRange(-1,0,100), LinRange(1,2,100)]
        t2 = (t[1:end-1]+t[2:end])/2
        y1 = eval_bspline.(1, t, T)
        y2 = eval_bspline_derivative.(1, 1, t2, T)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end

    t0 = LinRange(-1,7,100)
    for D in 0:6
        @test norm(eval_bspline_derivative.(D, D+1, t0, T))==0
    end

    t0 = -2.25:1/(1<<4):10
    t1 = (t0[1:end-1]+t0[2:end])/2
    t2 = (t1[1:end-1]+t1[2:end])/2
    t3 = (t2[1:end-1]+t2[2:end])/2
    t4 = (t3[1:end-1]+t3[2:end])/2
    t5 = (t4[1:end-1]+t4[2:end])/2

    for D in 3:6
        y0 = eval_bspline_derivative.(D, 0, t0, T)
        y2 = eval_bspline_derivative.(D, 2, t2, T)
        @test .5 > norm(diff(diff(y0)./diff(t0))./diff(t1) - y2)
    end

    for D in 4:6
        y0 = eval_bspline_derivative.(D, 0, t0, T)
        y2 = eval_bspline_derivative.(D, 3, t3, T)
        @test .5 > norm(  diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2) - y2)
    end

    for D in 5:6
        y0 = eval_bspline_derivative.(D, 0, t0, T)
        y2 = eval_bspline_derivative.(D, 4, t4, T)
        @test .5 > norm(  diff(diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2))./diff(t3) - y2)
    end

    for D in 6:8
        y0 = eval_bspline_derivative.(D, 0, t0, T)
        y2 = eval_bspline_derivative.(D, 5, t5, T)
        @test .5 > norm(  diff(diff(diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2))./diff(t3))./diff(t4) - y2)
    end

    for D in 1:6
        @test Int.(eval_bspline_derivative.(D, D, -.5:D+2, T)) == vcat(0,[(-1)^k*binomial(D,k) for k in 0:D],0)
    end
end

function allocation_test()
    T = Float64
    for f in (eval_bspline, eval_centered_bspline)
        for K in 0:10
            for x in LinRange(-10,10,30)
                @test @timed(f(K, x, T))[3] <= 32
            end
        end
    end
    for f in (eval_periodic_bspline, eval_periodic_centered_bspline)
        for K in 0:10
            for x in LinRange(-10,10,30)
                # Z = @timed(f(K, x, 1.0, T))[3]
                # @show K, Z
                @test @timed(f(K, x, 1.0, T))[3] <= 32
            end
        end
    end
    for K in 0:10, D in 1:K
        f = t -> eval_bspline_derivative(K, D, t, Float64)
        for x in LinRange(-10,10,2)
            f(x)
        end
        for x in LinRange(-10,10,30)
            @test @timed(f(x))[3] <= 32
        end
    end
end

@testset "B-spline implementation" begin
    for T in (Float64, BigFloat)
        @testset "elementary properties" begin
            elementary_props_of_splines_test(T)
        end
        @testset "periodic B splines"  begin
            periodic_bsplines_test(T)
        end
        @testset "centered B splines"  begin
            centered_bsplines_test(T)
        end
        @testset "derivative of B splines"  begin
            derivative_test(T)
        end
    end
    @testset "integrals of B splines"  begin
        test_spline_integrals()
    end
    @testset "allocation of evaluation of B splines"  begin
        allocation_test()
    end
end
