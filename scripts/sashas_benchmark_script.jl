using Revise, Groebner, Nemo, RationalUnivariateRepresentation, Statistics, TimerOutputs
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

coeff_size = rur -> Int(round(maximum(f -> maximum(c -> log2(abs(numerator(c))) + log2(abs(denominator(c))), f; init=-1), rur; init=-1)))

coeff_size2 = rur -> begin
    rur_z = map(f -> f*lcm(map(denominator, f)), rur)
    low=Int(round(minimum(f -> minimum(c -> log2(max(abs(numerator(c)),1)), f; init=-1), rur_z; init=-1)))
    med=Int(round(median(map(f -> length(f) == 0 ? -1 : median(map(c -> log2(max(abs(numerator(c)),1)), f)), rur_z))))
    up=Int(round(maximum(f -> maximum(c -> log2(max(abs(numerator(c)),1)), f; init=-1), rur_z; init=-1)))
    tot=Int(round(sum(f -> sum(c -> log2(max(abs(numerator(c)),1)), f; init=0), rur_z; init=0)))
    (low, med, up, tot)
end

function from_template(strat, rur, sep, time, m, to)
    if all(iszero, sep) || rur == nothing
        return """
	search_strategy=:$strat
	    FAILED"""
    end

    io = IOBuffer()
    show(io, to, compact = true, allocations = false, linechars = :ascii)
    
    """
search_strategy=:$strat
    sep. form    $sep   $(sum(!iszero, sep)) / $(length(sep))
    coeff size   $(digitsep(coeff_size(rur))) bits
    coeff size2  $(join(digitsep.(coeff_size2(rur)), ", ")) bits
    sparsity     $(round(sparsity(m), digits=2))
    time         $(round(time, digits=2)) s
    
$(String(take!(io)))

    """
end

for name in [
        "caprasse",
        "phuoc1",
        # "goodwin",
        # "crn", 	# shape
        # "seir36",
        # "crauste-1",	# shape
        # "crauste-2",
        # "fab_4", 
        # "noon6", 
        # "reimer6",
        # "robot",	# shape
        # "chandra10",	# shape
        # "chandra11",	# shape
        # "eco10",	# shape
        # "eco11",	# shape
        # "katsura11",
        # "Ch6_sq",
        # "Ka6_sq",
        # "No5_sq",
        # "Re5_sq",
        # "Ro5_sq",
        # "root7",
        # "crauste-3",	# too large
    ]

    include((@__DIR__) * "/../Data/Systems/$name.jl")

    infos = guess_infos(sys)

    # ring_zp, (x_zp..., T) = polynomial_ring(AbstractAlgebra.GF(2^30+3), vcat(map(string, gens(R)), "T"), internal_ordering=:degrevlex)
    ring_zp, x_zp = polynomial_ring(AbstractAlgebra.GF(2^30+3), map(string, gens(R)), internal_ordering=:degrevlex)
    sys_zp = map(f -> evaluate(map_coefficients(c -> base_ring(ring_zp)(BigInt(numerator(c))) // base_ring(ring_zp)(BigInt(denominator(c))), f), x_zp), sys)
    # sep1 = rand(1:100, length(vcat(x_zp, T)))
    # sys_zp = vcat(sys_zp, sum(sep1 .* vcat(x_zp, T)))
    matrices = multiplication_matrices(sys_zp)

    # reset_timer!(RationalUnivariateRepresentation.to)
    # time1 = @elapsed flag1, rur1 = RationalUnivariateRepresentation.rur_core(sys_zp)
    # to1 = deepcopy(RationalUnivariateRepresentation.to)
    
    reset_timer!(RationalUnivariateRepresentation.to)
    time1 = @elapsed rur1, sep1 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:current)
    to1 = deepcopy(RationalUnivariateRepresentation.to)

    reset_timer!(RationalUnivariateRepresentation.to)
    RationalUnivariateRepresentation._BOUND[] = 10
    sep2 = zeros(Int, length(matrices))
    time2, rur2 = nothing, nothing
    try
    	time2 = @elapsed rur2, sep2 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:random)
    catch e
	bound = RationalUnivariateRepresentation._BOUND[]
	@info "RUR failed with random linear form -$(bound)..$(bound) \\ 0"
    end
    to2 = deepcopy(RationalUnivariateRepresentation.to)

    reset_timer!(RationalUnivariateRepresentation.to)
    RationalUnivariateRepresentation._BOUND[] = 100
    sep3 = zeros(Int, length(matrices))
    time3, rur3 = nothing, nothing
    try
        time3 = @elapsed rur3, sep3 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:random)
    catch e
        bound = RationalUnivariateRepresentation._BOUND[]
        @info "RUR failed with random linear form -$(bound)..$(bound) \\ 0"
    end
    to3 = deepcopy(RationalUnivariateRepresentation.to)

    # reset_timer!(RationalUnivariateRepresentation.to)
    # time4 = @elapsed rur4, sep4 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:l0_norm)
    # to4 = deepcopy(RationalUnivariateRepresentation.to)

    # reset_timer!(RationalUnivariateRepresentation.to)
    # time5 = @elapsed rur5, sep5 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:mron_0l)
    # to5 = deepcopy(RationalUnivariateRepresentation.to)
    
    # m1 = zeros(length(sep1), length(sep1)) 
    # sum(sep1 .* matrices)
    m1 = sum(sep1 .* matrices)
    m2 = sum(sep2 .* matrices)
    m3 = sum(sep3 .* matrices)
    # m4 = sum(sep4 .* matrices)
    # m5 = sum(sep5 .* matrices)

    # results = """
    # =================================
    # $name
    #     n            $(length(gens(R)))
    #     K            $(base_ring(ring_zp))
    #     D            $(infos.quotient)
    #     d            $(infos.minpoly)
    #     type         $(infos.type)

    # $(from_template(:random100, rur1, sep1, time1, m1, to1))
    # """

    #         sparsity M_x $(map(x -> round(x, digits=2), map(sparsity, matrices)))

    results = """
    =================================
    $name
        n            $(length(gens(R)))
        D            $(infos.quotient)
        d            $(infos.minpoly)
        type         $(infos.type)
        sparsity M_x $(map(x -> round(x, digits=2), map(sparsity, matrices)))

    $(from_template(:current, rur1, sep1, time1, m1, to1))

    $(from_template(:random10, rur2, sep2, time2, m2, to2))

    $(from_template(:random100, rur3, sep3, time3, m3, to3))
        """
    #  

    # $(from_template(:l0_norm, rur4, sep4, time4, m4, to4))

    # $(from_template(:mron_0l, rur5, sep5, time5, m5, to5))
    println(results)
    open((@__DIR__) * "/results.txt", "a") do file println(file, results) end
end

