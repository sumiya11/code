using AbstractAlgebra, RationalUnivariateRepresentation

include((@__DIR__) * "/../Data/Systems/Ro5_sq.jl")

rur1, sep1 = rur(sys, get_separating_element=true, search_strategy=:current);

rur2, sep2 = rur(sys, get_separating_element=true, search_strategy=:deterministic);

rur3, sep3 = rur(sys, get_separating_element=true, search_strategy=:random);

