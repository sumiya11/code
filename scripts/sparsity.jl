using Groebner, Nemo, RationalUnivariateRepresentation
import AbstractAlgebra

function multiplication_matrices(J)
    J = Groebner.groebner(J)
    basis = sort(Groebner.quotient_basis(J))
    dim = length(basis)
    S = Nemo.matrix_space(Nemo.base_ring(J[1]), dim, dim)
    matrices = []
    @info "Dim is $dim"
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

    ring_zp, x_zp = polynomial_ring(GF(2^30+3), map(string, gens(R)), internal_ordering=:degrevlex)
    sys_zp = map(f -> evaluate(map_coefficients(c -> base_ring(ring_zp)(c), f), x_zp), sys)
    matrices = multiplication_matrices(sys_zp)
    
    rur1, sep1 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:current)
    rur2, sep2 = zdim_parameterization(sys, get_separating_element=true, search_strategy=:random)
    
    m1 = sum(sep1 .* matrices)
    m2 = sum(sep2 .* matrices)
    
    println("""
        =================================
        $name
        
        Multiplication matrices sparsity:
        $(map(x -> round(x, digits=3), map(sparsity, matrices)))
    
        search_strategy=:current
            sep. form    $sep1
            sparsity     $(round(sparsity(m1), digits=3))
    
        search_strategy=:random
            sep. form    $sep2
            sparsity     $(round(sparsity(m2), digits=3))
        """)
end

#=
=================================
Ro5_sq

Multiplication matrices sparsity:
[0.062, 0.043, 0.031, 0.023, 0.014]

search_strategy=:current
    sep. form    [2, -2, 0, -1, 1]
    sparsity     0.068

search_strategy=:random
    sep. form    [26, -10, -5, 41, 35]
    sparsity     0.068

=================================
caprasse

Multiplication matrices sparsity:
[0.094, 0.087, 0.083, 0.073]

search_strategy=:current
    sep. form    [2, 0, -1, 1]
    sparsity     0.196

search_strategy=:random
    sep. form    [26, -10, -5, 41]
    sparsity     0.24

=================================
fab_4

Multiplication matrices sparsity:
[0.23, 0.227, 0.209, 0.102]

search_strategy=:current
    sep. form    [-1, -1, -1, 1]
    sparsity     0.236

search_strategy=:random
    sep. form    [26, -10, -5, 41]
    sparsity     0.236

=================================
noon6

Multiplication matrices sparsity:
[0.035, 0.03, 0.029, 0.027, 0.024, 0.02]

search_strategy=:current
    sep. form    [-11, 5, -3, 0, -1, 1]
    sparsity     0.09

search_strategy=:random
    sep. form    [26, -10, -5, 41, 35, -67]
    sparsity     0.099

=================================
reimer6

Multiplication matrices sparsity:
[0.353, 0.357, 0.269, 0.244, 0.175, 0.142]

search_strategy=:current
    sep. form    [1, -1, 0, 0, -1, 1]
    sparsity     0.426

search_strategy=:random
    sep. form    [26, -10, -5, 41, 35, -67]
    sparsity     0.434
=#
