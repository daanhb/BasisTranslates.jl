
# simple implementation of the AZ algorithm for fixed rank
function az(A, Zstar, b; rank_estimate::Int = 30)
    m, n = size(A)
    R = rand(eltype(b), n, rank_estimate)
    Col = zeros(eltype(b), m, rank_estimate)
    for i in 1:rank_estimate
        r = R[:,i]
        u = real(A*r - A*(Zstar*(A*r)))
        Col[:,i] = u
    end
    u1 = Col \ real(b - A*(Zstar*b))
    x1 = R*u1
    x2 = real(Zstar*(b - A*x1))
    x = x1+x2
    x
end
