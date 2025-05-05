    # ring_zp, (x_zp..., T) = polynomial_ring(AbstractAlgebra.GF(2^30+3), vcat(map(string, gens(R)), "T"), internal_ordering=:degrevlex)
    ring_zp, x_zp = polynomial_ring(AbstractAlgebra.GF(2^30+3), map(string, gens(R)), internal_ordering=:degrevlex)
    sys_zp = map(f -> evaluate(map_coefficients(c -> base_ring(ring_zp)(BigInt(numerator(c))) // base_ring(ring_zp)(BigInt(denominator(c))), f), x_zp), sys)
    # sep1 = rand(1:100, length(vcat(x_zp, T)))
    # sys_zp = vcat(sys_zp, sum(sep1 .* vcat(x_zp, T)))
    matrices = multiplication_matrices(sys_zp)

    # reset_timer!(RationalUnivariateRepresentation.to)
    # time1 = @elapsed flag1, rur1 = RationalUnivariateRepresentation.rur_core(sys_zp)
    # to1 = deepcopy(RationalUnivariateRepresentation.to)
    

# results = """
    # =================================
    # $name
    #     n            $(length(gens(R)))
    #     K            $(base_ring(ring_zp))
    #     D            $(infos.quotient)
    #     d            $(infos.minpoly)
    #     type         $(infos.type)

    # $(from_template(:random100, rur1, sep1, time1, m1, to1))
    # ""
    #
    # "
