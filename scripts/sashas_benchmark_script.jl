using AbstractAlgebra, RationalUnivariateRepresentation

for name in [
        # non-trivial separating form
        "Ro5_sq",
        "caprasse", 
        "fab_4", 
        "noon6", 
        "reimer6",
        # "root7"   fails for search_strategy=:random
    ]

    include((@__DIR__) * "/../Data/Systems/$name.jl")
    
    time1 = @elapsed rur1, sep1 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:current)
    
    time2 = @elapsed rur2, sep2 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:random)

    coeff_size = rur -> Int(round(maximum(f -> maximum(c -> log2(abs(numerator(c))) + log2(abs(denominator(c))), f), rur)))
    println("""
    =================================
    $name

    search_strategy=:current
        sep. form   $sep1
        coeff size  $(coeff_size(rur1)) bits
        time        $(round(time1, digits=2)) s

    search_strategy=:random
        sep. form   $sep2
        coeff size  $(coeff_size(rur2)) bits
        time        $(round(time2, digits=2)) s
        """)
end

#=
=================================
caprasse

search_strategy=:current
    sep. form   [2, 0, -1, 1]
    coeff size  56 bits
    time        0.01 s

search_strategy=:random
    sep. form   [26, -10, -5, 41]
    coeff size  177 bits
    time        0.01 s

=================================
fab_4

search_strategy=:current
    sep. form   [-1, -1, -1, 1]
    coeff size  28056 bits
    time        175.3 s

search_strategy=:random
    sep. form   [26, -10, -5, 41]
    coeff size  32837 bits
    time        170.99 s

=================================
noon6

search_strategy=:current
    sep. form   [-11, 5, -3, 0, -1, 1]
    coeff size  4087 bits
    time        7.88 s

search_strategy=:random
    sep. form   [26, -10, -5, 41, 35, -67]
    coeff size  6034 bits
    time        7.67 s

=================================
reimer6

search_strategy=:current
    sep. form   [1, -1, 0, 0, -1, 1]
    coeff size  1924 bits
    time        1.89 s

search_strategy=:random
    sep. form   [26, -10, -5, 41, 35, -67]
    coeff size  5638 bits
    time        4.25 s
=#

