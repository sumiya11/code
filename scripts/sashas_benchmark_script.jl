using Groebner, Nemo, RationalUnivariateRepresentation
import AbstractAlgebra

function multiplication_matrices(J)
    J = Groebner.groebner(J)
    basis = sort(Groebner.quotient_basis(J))
    dim = length(basis)
    S = Nemo.matrix_space(Nemo.base_ring(J[1]), dim, dim)
    matrices = []
    println("Dim is $dim")
    for v in gens(parent(first(J)))
        print(v, ",")
        M = zero(S)
        images = Groebner.normalform(J, v .* basis)
        for i in 1:length(basis)
            for (j, base_vec) in enumerate(basis)
                M[i, j] = coeff(images[i], base_vec)
            end
        end
        push!(matrices, M)
    end
    println()
    matrices
end

sparsity = m -> sum(!iszero, m) / prod(size(m))
coeff_size = rur -> Int(round(maximum(f -> maximum(c -> log2(abs(numerator(c))) + log2(abs(denominator(c))), f), rur)))

for name in [
        # "caprasse", 
        # "fab_4", 
        "noon6", 
        "reimer6",
        # "robot",
        # "chandra10",
        # "chandra11",
        # "eco10",
        # "eco11",
        # "katsura11",
        # "Ch6_sq",
        # "Ka6_sq",
        # "No5_sq",
        # "Re5_sq",
        "Ro5_sq",
        # "root7",   fails for search_strategy=:random
    ]

    include((@__DIR__) * "/../Data/Systems/$name.jl")

    infos = guess_infos(sys)

    ring_zp, x_zp = polynomial_ring(GF(2^30+3), map(string, gens(R)), internal_ordering=:degrevlex)
    sys_zp = map(f -> evaluate(map_coefficients(c -> base_ring(ring_zp)(c), f), x_zp), sys)
    matrices = multiplication_matrices(sys_zp)
    
    time1 = @elapsed rur1, sep1 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:current)
    
    time2 = @elapsed rur2, sep2 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:random)

    time3 = @elapsed rur3, sep3 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:l0_norm)

    m1 = sum(sep1 .* matrices)
    m2 = sum(sep2 .* matrices)
    m3 = sum(sep3 .* matrices)
    
    results = """
    =================================
    $name
       D            $(infos.quotient)
       d            $(infos.minpoly)
       sparsity     $(map(x -> round(x, digits=2), map(sparsity, matrices)))
       type         $(infos.type)

    search_strategy=:current
        sep. form   $sep1
        coeff size  $(coeff_size(rur1)) bits
        sparsity    $(round(sparsity(m1), digits=2))
        time        $(round(time1, digits=2)) s

    search_strategy=:random
        sep. form   $sep2
        coeff size  $(coeff_size(rur2)) bits
        sparsity    $(round(sparsity(m2), digits=2))
        time        $(round(time2, digits=2)) s

    search_strategy=:l0_norm
        sep. form   $sep3
        coeff size  $(coeff_size(rur3)) bits
        sparsity    $(round(sparsity(m3), digits=2))
        time        $(round(time3, digits=2)) s
        """
    println(results)
    open((@__DIR__) * "/results.txt", "a") do file println(file, results) end
end

