
# starting Julia with "julia --threads=N"

using Pkg; Pkg.develop(path="../Julia/RationalUnivariateRepresentation.jl/")
Pkg.status()

using InteractiveUtils, Base.Threads
using Groebner, RationalUnivariateRepresentation, ThreadPinning

# For better reproducibility
pinthreads(:cores)

versioninfo()
#threadinfo()

for name in ["schwarz11"] # ["caprasse", "eco10", "chandran9", "noon6", "reimer6"]
    include("../Data/Systems/$name.jl")
    @info "name = $name"
    @time rur1 = zdim_parameterization(sys, verbose=false, threads=1)
    @time rur2 = zdim_parameterization(sys, verbose=false, threads=4)
    @time rur3 = zdim_parameterization(sys, verbose=false, threads=8)
    @time rur4 = zdim_parameterization(sys, verbose=false, threads=16)
    @assert rur1 == rur2 == rur3 == rur4 
end

#=
threads=N

kat 11
    total:
        N=1          358 s
        N=4          121 s
        N=8           97 s
        N=16          87 s

    prepare_system:    7 s
    learn_zdim_quo:   12 s

eco 12
    total:
        N=1          610 s
        N=4          221 s
        N=8          173 s
        N=16         130 s

    prepare_system:   10 s
    learn_zdim_quo:   23 s

schwarz11
    total:
        N=1         4221 s
        N=4         1837 s
        N=8          552 s
        N=16         413 s

    prepare_system:  110 s
    learn_zdim_quo:   28 s

chandran9
    total:
        N=1           70 s
        N=4           30 s
        N=8           29 s
        N=16          22 s
        
    prepare_system:    0 s
    learn_zdim_quo:    1 s
=#
