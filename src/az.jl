
function tsvd(A, b; tol)
    u,s,v = svd(A; full=false)
    I = findlast(s/s[1] .>= tol)
    if I < size(A,2)
        u = u[:,1:I]
        v = ((v')[1:I,:])'
        s = s[1:I]
    end
    sinv = s.^(-1)
    v*(sinv .* (u'*b))
end

# simple implementation of the AZ algorithm for fixed rank
function az(A, Zstar, b; rank_estimate::Int = 30, tol = 1e-12)
    m, n = size(A)
    R = rand(eltype(b), n, rank_estimate)
    Col = zeros(eltype(b), m, rank_estimate)
    for i in 1:rank_estimate
        r = R[:,i]
        u = real(A*r - A*(Zstar*(A*r)))
        Col[:,i] = u
    end
    u1 = tsvd(Col, real(b - A*(Zstar*b)); tol)
    x1 = R*u1
    x2 = real(Zstar*(b - A*x1))
    x = x1+x2
    x
end
