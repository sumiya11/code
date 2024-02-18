# starting Julia with "julia --threads=N --procs=(N-1)"
using Distributed

@everywhere using Pkg; Pkg.develop(path=(@__DIR__)*"/../Julia/RationalUnivariateRepresentation.jl/")
@everywhere using RationalUnivariateRepresentation

Pkg.status()

using InteractiveUtils, Base.Threads
using Groebner, ThreadPinning

versioninfo()

for name in ["caprasse", "eco10", "chandran9", "noon6", "reimer6"]
    include("../Data/Systems/$name.jl")
    @info "name = $name"
    @time rur1 = zdim_parameterization(sys, parallelism=:serial)
    @time rur2 = zdim_parameterization(sys, parallelism=:multithreading)
    @time rur3 = zdim_parameterization(sys, parallelism=:multiprocessing, procs=Distributed.procs())
    @assert rur1 == rur2 == rur3 
end

#=

=#
