function az_approximate2(basis::UnitPeriodicTranslates, f, t1, t2, osf;
        options...)

    n = length(basis)
    tt = oversampled_grid(basis; osf)

    # Compute the restriction R from the large grid to the subgrid
    I1 = findfirst(tt .>= t1)
    I2 = findlast(tt .<= t2)
    R = RestrictionArray(length(tt), I1:I2)
    E = R'
    tt_int = R*tt

    b = f.(tt_int)

    # Compute the factorization of A using its block-circulant structure
    Acirc = matrix(evaluation(basis, tt))
    F = factorize(Acirc)

    # Acirc2 = ArrayOperator(F.Π')*ArrayOperator(F.P)*ArrayOperator(F.D)*ArrayOperator(F.F')
    # Ae = ArrayOperator(E)*ArrayOperator(R)*Acirc2
    Ae = ArrayOperator(E)*ArrayOperator(R)*ArrayOperator(Acirc)
    Be = ArrayOperator(F.P) * Ae * ArrayOperator(F.F)
    be = F.P * (E*b)

    AZ_A = Be
    AZ_Zstar = ArrayOperator(F.Dpinv)
    y = az(AZ_A, AZ_Zstar, be; rank_estimate=90, options...)

    c = F.F * y
    Expansion(basis, c)
end

function az_approximate3(basis::UnitPeriodicTranslates, f, t1, t2, osf;
        options...)

    n = length(basis)
    tt = oversampled_grid(basis; osf)

    # Compute the restriction R from the large grid to the subgrid
    I1 = findfirst(tt .>= t1)
    I2 = findlast(tt .<= t2)
    R = RestrictionArray(length(tt), I1:I2)
    E = R'
    tt_int = R*tt

    b = f.(tt_int)

    # Compute the factorization of A using its block-circulant structure
    Acirc = matrix(evaluation(basis, tt))
    F = factorize(Acirc)

    Ae = ArrayOperator(E)*ArrayOperator(R)*ArrayOperator(Acirc)
    Be = ArrayOperator(F.P) * Ae * ArrayOperator(F.F)
    be = F.P * (E*b)

    AZ_A = Be
    AZ_Zstar = ArrayOperator(F.Dpinv)
    # y = az(AZ_A, AZ_Zstar, be; rank_estimate=60, options...)
    AZ_A0 = ArrayOperator(R)*ArrayOperator(Acirc) * ArrayOperator(F.F)
    AZ_Z = ArrayOperator(F.Dpinv)
    y = AZalgorithm.az0(AZ_A, AZ_Zstar, b, AZ_A0, ArrayOperator(F.P), ArrayOperator(E); rank_estimate=60, options...)

    c = F.F * y
    Expansion(basis, c)
end
