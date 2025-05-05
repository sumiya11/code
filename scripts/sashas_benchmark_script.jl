using Revise, Groebner, Nemo, RationalUnivariateRepresentation, Statistics, TimerOutputs, Dates, JLD2, PrettyTables
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

bitsize_zz(c) = (@assert denominator(c) == 1; log2(abs(max(1, numerator(c)))))
bitsize_qq(c) = log2(abs(max(1, numerator(c)))) + log2(abs(denominator(c)))

coeff_size_qq = rur -> begin
    low=Int(round(minimum(f -> minimum(bitsize_qq, f; init=-1), rur; init=-1)))
    med=Int(round(median(map(f -> length(f) == 0 ? -1 : median(map(bitsize_qq, f)), rur))))
    up=Int(round(maximum(f -> maximum(bitsize_qq, f; init=-1), rur; init=-1)))
    tot=Int(round(sum(f -> sum(bitsize_qq, f; init=0), rur; init=0)))
    (min=low, med=med, max=up, tot=tot)
end

coeff_size_zz = rur -> begin
    rur_z = map(f -> f*lcm(map(denominator, f)), rur)
    low=Int(round(minimum(f -> minimum(bitsize_zz, f; init=-1), rur_z; init=-1)))
    med=Int(round(median(map(f -> length(f) == 0 ? -1 : median(map(bitsize_zz, f)), rur_z))))
    up=Int(round(maximum(f -> maximum(bitsize_zz, f; init=-1), rur_z; init=-1)))
    tot=Int(round(sum(f -> sum(bitsize_zz, f; init=0), rur_z; init=0)))
    (min=low, med=med, max=up, tot=tot)
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
    sep. form     $sep   $(sum(!iszero, sep)) / $(length(sep))
    coeff sz qq   $(join(string.(keys(coeff_size_qq(rur))) .* "=" .* digitsep.(values(coeff_size_qq(rur))), ", ")) bits
    coeff sz zz   $(join(string.(keys(coeff_size_zz(rur))) .* "=" .* digitsep.(values(coeff_size_zz(rur))), ", ")) bits
    sparsity M_t  $(round(sparsity(m), digits=2))
    time          $(round(time, digits=2)) s
    
$(String(take!(io)))

    """
end

function main()
mkpath(joinpath(@__DIR__, "results"))
ID = Dates.format(now(), "yyyy_mm_dd") * "-" * string(length(filter(file -> endswith(file, r"-[0-9]+"), readdir(joinpath(@__DIR__, "results")))))
mkpath(joinpath(@__DIR__, "results", ID))

@info "Writing results to $(joinpath(@__DIR__, "results", ID))"

table = []

for name in [
        # ============= From files =============
        "caprasse",
        # "phuoc1",
        "goodwin",
        "crn", 	        # shape
        "seir36",
        "crauste-1",	# shape
        "crauste-2",
        "schwarz11",
        # "crauste-3",	# too large
        "fab_4", 
        "robot",	    # shape
        "Ch6_sq",
        "Ka6_sq",
        "No5_sq",
        "Re5_sq",
        "Ro5_sq",
        # ============= From Groebner =============
        "Noon_5" => Groebner.Examples.noonn(5),
        "Noon_6" => Groebner.Examples.noonn(6),
        "Noon_7" => Groebner.Examples.noonn(7),
        "Chandra_9" => Groebner.Examples.chandran(9),
        "Chandra_10" => Groebner.Examples.chandran(10),
        "Chandra_11" => Groebner.Examples.chandran(11),
        "Katsura_10" => Groebner.Examples.katsuran(9),
        "Katsura_11" => Groebner.Examples.katsuran(10),
        "Katsura_12" => Groebner.Examples.katsuran(11),
        "Reimer_6" => Groebner.Examples.reimern(6),
        "Reimer_7" => Groebner.Examples.reimern(7),
        "Eco_11" => Groebner.Examples.econ(11),
        "Eco_12" => Groebner.Examples.econ(12),
        "Eco_13" => Groebner.Examples.econ(13),
        "Henrion_6" => Groebner.Examples.henrion6(),
    ]

    if !(name isa Pair)
        sys = include(joinpath(@__DIR__, "../Data/Systems/$name.jl"))
    else
        name, sys = name
    end
    
    infos = guess_infos(sys)

    ring_zp, x_zp = polynomial_ring(AbstractAlgebra.GF(2^30+3), map(string, gens(parent(sys[1]))), internal_ordering=:degrevlex)
    sys_zp = map(f -> evaluate(map_coefficients(c -> base_ring(ring_zp)(BigInt(numerator(c))) // base_ring(ring_zp)(BigInt(denominator(c))), f), x_zp), sys)
    
    matrices = multiplication_matrices(sys_zp)

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

    reset_timer!(RationalUnivariateRepresentation.to)
    time4 = @elapsed rur4, sep4 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:mron_0l)
    to4 = deepcopy(RationalUnivariateRepresentation.to)

    m1 = sum(sep1 .* matrices)
    m2 = sum(sep2 .* matrices)
    m3 = sum(sep3 .* matrices)
    m4 = sum(sep4 .* matrices)
    # m5 = sum(sep5 .* matrices)

    results = """
    =================================
    $name
        n             $(length(gens(parent(sys[1]))))
        K             $(base_ring(parent(sys[1])))
        D             $(infos.quotient)
        d             $(infos.minpoly)
        type          $(infos.type)
        sparsity M_x  $(map(x -> round(x, digits=2), map(sparsity, matrices)))

    $(from_template(:current, rur1, sep1, time1, m1, to1))

    $(from_template(:random10, rur2, sep2, time2, m2, to2))

    $(from_template(:random100, rur3, sep3, time3, m3, to3))

    $(from_template(:mron_0l, rur4, sep4, time4, m4, to4))
        """
    # $(from_template(:l0_norm, rur4, sep4, time4, m4, to4))
    # $(from_template(:mron_0l, rur5, sep5, time5, m5, to5))
    
    println(results)
    mkpath(joinpath(@__DIR__, "results", ID, name))
    open(joinpath(@__DIR__, "results", ID, name, "results.txt"), "w") do file println(file, results) end
    for (strategy, rur, sep, time) in [(:current, rur1, sep1, time1), (:random10, rur2, sep2, time2), (:random100, rur3, sep3, time3), (:mron_0l, rur4, sep4, time4)]
        jldsave(joinpath(@__DIR__, "results", ID, name, "rur_$strategy.jld2"); rur=rur, sep=sep)
        push!(table, [name, infos.quotient, infos.minpoly, infos.type, strategy, time])
    end
end

pretty_table(permutedims(reduce(hcat, table)), title="Summary", header=["Name", "D", "d", "Type", "Strategy", "Time"], limit_printing=false)
end

main()

