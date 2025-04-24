using Groebner, Nemo, RationalUnivariateRepresentation, Statistics
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

# from https://github.com/IainNZ/Humanize.jl/blob/master/README.md
function digitsep(value::Integer; seperator="_", per_separator=3)
    isnegative = value < zero(value)
    value = string(abs(value))  # Stringify, no seperators.
    # Figure out last character index of each group of digits.
    group_ends = reverse(collect(length(value):-per_separator:1))
    groups = [value[max(end_index - per_separator + 1, 1):end_index]
              for end_index in group_ends]
    return (isnegative ? "-" : "") * join(groups, seperator)
end

sparsity = m -> sum(!iszero, m) / prod(size(m))

coeff_size = rur -> Int(round(maximum(f -> maximum(c -> log2(abs(numerator(c))) + log2(abs(denominator(c))), f), rur)))

coeff_size2 = rur -> begin
    rur_z = map(f -> f*lcm(map(denominator, f)), rur)
    low=Int(round(minimum(f -> minimum(c -> log2(max(abs(numerator(c)),1)), f), rur_z)))
    med=Int(round(median(map(f -> median(map(c -> log2(abs(numerator(c))), f)), rur_z))))
    up=Int(round(maximum(f -> maximum(c -> log2(abs(numerator(c))), f), rur_z)))
    tot=Int(round(sum(f -> sum(c -> log2(max(abs(numerator(c)),1)), f), rur_z)))
    (low, med, up, tot)
end

function from_template(strat, rur, sep, time, m)
    """
search_strategy=:$strat
    sep. form    $sep
    coeff size   $(digitsep(coeff_size(rur))) bits
    coeff size2  $(join(digitsep.(coeff_size2(rur)), ", ")) bits
    sparsity     $(round(sparsity(m), digits=2))
    time         $(round(time, digits=2)) s
    """
end

for name in [
        # "caprasse",
        # "goodwin",
        # "crn",
        # "seir36",
        "crauste-1",
        "crauste-2",
        # "crauste-3",
        "fab_4", 

        
        "noon6", 
        "reimer6",
        "robot",
        "chandra10",
        "chandra11",
        "eco10",
        "eco11",
        "katsura11",
        "Ch6_sq",
        "Ka6_sq",
        # "No5_sq",
        "Re5_sq",
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

    time4 = @elapsed rur4, sep4 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:mron_0l)

    m1 = sum(sep1 .* matrices)
    m2 = sum(sep2 .* matrices)
    m3 = sum(sep3 .* matrices)
    m4 = sum(sep4 .* matrices)
    
    results = """
    =================================
    $name
        n            $(length(gens(R)))
        D            $(infos.quotient)
        d            $(infos.minpoly)
        type         $(infos.type)
        sparsity     $(map(x -> round(x, digits=2), map(sparsity, matrices)))

    $(from_template(:current, rur1, sep1, time1, m1))

    $(from_template(:random, rur2, sep2, time2, m2))

    $(from_template(:l0_norm, rur3, sep3, time3, m3))

    $(from_template(:mron_0l, rur4, sep4, time4, m4))
        """
    println(results)
    open((@__DIR__) * "/results.txt", "a") do file println(file, results) end
end

